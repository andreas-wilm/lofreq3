## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License

{.compile: "viterbi.c".}
{.passL: "-lm"}

# standard
import strformat
import strutils
import tables
import sequtils
import algorithm

# third party
import hts

# project specific
#/

# void
proc left_align_indels(sref: cstring, squery: cstring, slen: cint,
  new_state_seq: cstring) {.cdecl, importc: "left_align_indels".}


# returns shift
proc viterbi_c(sref: cstring, squery: cstring, bqual: ptr uint8,
  saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}


proc skipRead(rec: Record, skipSecondary: bool): bool =

  # skip secondary reads unless requested
  if rec.flag.secondary and skipSecondary:
    return true

  # skip unmapped reads
  if rec.flag.unmapped:
    return true

  # skip reads without indels
  var hasIndels = false
  for c in rec.cigar:
    if c.op in @[CigarOp.insert, CigarOp.deletion]:
      hasIndels = true
      break
  if not hasIndels:
    return true

  return false


# get reference context for aligned read
proc getRefContext(rec: Record, refSq: string, refPadding: int): string =
  var s = rec.start - refPadding
  if s < 0:
    s = 0
  var e = rec.stop + refPadding
  if e >= len(refSq):
    e = len(refSq)
  result = toUpperAscii(refSq[s..e])#toupper to avoid ref masking


# from https://github.com/brentp/bamject/blob/master/src/cigar.nim

const BAM_CIGAR_SHIFT = 4'u32
const BAM_CIGAR_STR = "MIDNSHP=XB"

var cigar_tab : array[128, int]
for i in 0..<128:
    cigar_tab[i] = -1
for i, c in BAM_CIGAR_STR:
    cigar_tab[c.int] = i


proc tocigar(cs:string): seq[CigarElement] =
  if cs.len > 4:
    result = newSeqOfCap[CigarElement](2)
  var off = 0
  while off < cs.len:
    var i = 0
    while cs[off + i].isdigit:
      i += 1
    var num = parseInt(cs[off..<off+i])
    var ops = cs[off+i]
    off += i + 1
    var el:uint32 = num.uint32 shl BAM_CIGAR_SHIFT
    if cigar_tab[ops.int] == -1:
      quit "unknown cigar op from " & cs & ": " & $ops
    el = el or cigar_tab[ops.int].uint32
    result.add(cast[CigarElement](el))

# end from

# fold unrolled cigarstring, e.g. MMMMDD to 4M2D
proc foldCigar(unfoldedCigar: string): string =
  var lastop = unfoldedCigar[0]
  var oplen = 1
  for i in countup(1, len(unfoldedCigar)-1):
    let op = unfoldedCigar[i]
    if op == lastop:
      oplen += 1
    else:
      result.add($oplen & lastop)
      oplen = 1
      lastop = op
  result.add($oplen & lastop)


proc createRealnRec(rec: Record, realnStart: int64, fullRealnCigar: Cigar): string =
    # We can't set BAM values in htsnim (argh), so convert to string / SAM
    # as in https://github.com/brentp/bamject/blob/master/src/bamject.nim
    var recSplit = rec.tostring().split('\t')
    recSplit[3] = $(realnStart + 1)# 1-based
    recSplit[5] = $fullRealnCigar

    # remove tags invalid after realignment.
    # FIXME ideally we should recompute them and set OA.
    # hope that fixmate fixes all mate references e.g. MC
    let tagsInvalidAfterRealn = @["MD", "AS", "NM"]
    var delIndices: seq[int]
    for i in countup(11, len(recSplit)-1):
      let fieldSplit = recSplit[i].split(":")
      if fieldSplit[0] in tagsInvalidAfterRealn:
        delIndices.add(i)
    for i, j in delIndices.pairs:
      recSplit.delete(j-i)
    result = recSplit.join("\t")


proc findSkipOps(rec: Record): (seq[CigarElement], int, seq[CigarElement], int) = 
    # find leading and trailing skip ops
    var leadingSkipOps, trailingSkipOps: seq[CigarElement]
    let skipOps = @[CigarOp.soft_clip, CigarOp.hard_clip, CigarOp.pad]
    for i in countup(0, len(rec.cigar)-1):
      var c = rec.cigar[i]
      if c.op notin skipOps:
        break
      leadingSkipOps.add(c)
    for i in countdown(len(rec.cigar)-1, len(leadingSkipOps)+1):
      var c = rec.cigar[i]
      if c.op notin skipOps:
        break
      trailingSkipOps.add(c)
    #stderr.writeLine("DEBUG cigar=" & $rec.cigar & " leadingSkipOps=" & $leadingSkipOps & " trailingSkipOps=" & $trailingSkipOps)

    # FIXME is first or last non skip pos is I/D warn, print and continue

    var leadingSoftClipLen = 0
    for c in leadingSkipOps:
      if c.op == CigarOp.soft_clip:
        leadingSoftClipLen.inc(c.len)
    var trailingSoftClipLen = 0
    for c in trailingSkipOps:
      # starting from lower index, but order doesn't matter here, since it's all skips
      if c.op == CigarOp.soft_clip:
        trailingSoftClipLen.inc(c.len)
    #stderr.writeLine("DEBUG leadingSoftClipLen=" & $leadingSoftClipLen & " trailingSoftClipLen=" & $trailingSoftClipLen)

    return (leadingSkipOps, leadingSoftClipLen, 
      trailingSkipOps, trailingSoftClipLen)


proc median(xs: seq[uint8]): uint8 =
  if len(xs) == 1: return xs[0]
  var ys = xs
  sort(ys, system.cmp[uint8])
  if len(ys) mod 2 == 0:
    # even number: return mean of the two elements in the middle
    result = (ys[len(ys) div 2] + ys[len(ys) div 2 - 1]) div 2;
  else:
    # odd number: return element in middle
     result = ys[len(ys) div 2];


proc medianQual(quals: seq[uint8]): uint8 =
  var nonQ2quals: seq[uint8]
  for q in quals:
    if q != 2:
      nonQ2quals.add(q)
  if len(nonQ2quals) == 0:
    return 0
  else:
    return median(nonQ2quals)


proc viterbi*(faFname: string, bamInFname: string, skipSecondary = true, refPadding = 10) =
  var fai: Fai
  # keeping all observerd reference sequences in memory for speedup
  var refs = initTable[string, string]()
  var iBam: Bam
  #var oBam: Bam

  if not open(fai, faFname):
    quit("Could not open reference sequence file: " & faFname)

  open(iBam, bamInFname, fai=faFname)
  
  stdout.write($iBam.hdr)
  #assert outFormat in @["BAM", "SAM", "CRAM"]
  #open(oBam, "-", mode="w" & outFormat, fai=faFname)
  #obam.write_header(iBam.hdr)
  
  for rec in iBam:
    var chrom = rec.chrom

    if skipRead(rec, skipSecondary):
      #obam.write(rec)
      echo $rec.tostring()
      continue

    # load reference if not cached
    if not refs.hasKey(chrom):
      #stderr.writeLine("DEBUG Loading " & chrom)
      refs[chrom] = fai.get(chrom)

    #stderr.writeLine("DEBUG incoming read = " & $rec.tostring())

    let refContext = getRefContext(rec, refs[chrom], refPadding)

    let (leadingSkipOps, leadingSoftClipLen, 
      trailingSkipOps, trailingSoftClipLen) = findSkipOps(rec)

    var query: string
    discard rec.sequence(query)
    let queryWOSoftClip = query[leadingSoftClipLen .. len(query)-1-trailingSoftClipLen]

    var bqual: seq[uint8]
    discard rec.base_qualities(bqual)
    var bqualWOSoftCLip = bqual[leadingSoftClipLen .. len(bqual)-1-trailingSoftClipLen]

    #echo "DEBUG: origQuery=" & query
    #echo "DEBUG: queryWOSoftClip=" & queryWOSoftClip# & " bqualWOSoftCLip=" & $bqualWOSoftCLip

    let q2def = cint(medianQual(bqualWOSoftCLip))
    # don't realign Q2 only
    if q2def == 0:
      echo $rec.tostring()
      continue

    let realnCigarRawOveralloc = newString(max(len(queryWOSoftClip), len(refContext)))
    # WARNING with newStringOfCap we get a segfault, but now newCigarRaw is oversized
    # printing works, but len is incorrect. creating a copy to a second string with $ doesn't help
    # FIXME what to do with shift?
    let shift = viterbi_c(refContext, queryWOSoftClip,
      cast[ptr uint8](addr(bqualWOSoftCLip[0])), realnCigarRawOveralloc, q2def)
    #stderr.writeLine("DEBUG realnCigarRawOveralloc ='" & $realnCigarRawOveralloc & "'")

    # realnCigarRaw is overallocated and len doesn't work!?
    # so let's fix this now
    var endIdx = len(realnCigarRawOveralloc)-1# this is only needed because the rawCigar is overallocated and len doesn't work!?
    for i, opChar in realnCigarRawOveralloc.pairs:
      if not opChar.isUpperAscii:
        endIdx = i-1
        break
    var realnCigarRaw = realnCigarRawOveralloc[0..endIdx]

    let realnCigar = toCigar(foldCigar(realnCigarRaw))
    #stderr.writeLine("DEBUG realnCigar=" & $realnCigarRaw)
    var fullRealnCigarSeq = concat(leadingSkipOps, realnCigar, trailingSkipOps)
    # https://github.com/brentp/hts-nim/commit/0bf2683b17f5ccc27a7af4731f187452b287cb61
    GC_ref(fullRealnCigarSeq)
    let fullRealnCigar = newCigar(fullRealnCigarSeq)
    #stderr.writeLine("DEBUG fullRealnCigar=" & $fullRealnCigar)

    # check if start pos was shifted
    var realnStart = rec.start
    var lower = rec.start - refPadding;
    if lower < 0:
      lower = 0
    if shift - (rec.start - lower) != 0:
      realnStart = rec.start + (shift - (rec.start - lower))
      #stderr.writeLine("DEBUG: new realnStart = " & $realnStart & " lower=" & $lower & " shift=" & $shift)
    #else:
      #stderr.writeLine("DEBUG: new=old realnStart = " & $realnStart)
      
    echo createRealnRec(rec, realnStart, fullRealnCigar)

    #oBam.close()

    GC_unref(fullRealnCigarSeq)

  stderr.writeLine("WARNING: MC tag in realigned mates invalid")


when isMainModule:
  block:
    var quals: seq[uint8]
    quals = @[2u8, 2u8, 2u8, 2u8]
    doAssert medianQual(quals) == 0
    quals = @[2u8, 10u8, 20u8, 40u8]
    doAssert medianQual(quals) == 20
    quals = @[2u8, 10u8, 20u8, 40u8, 100u8]
    doAssert medianQual(quals) == 30

  block:  
    doAssert foldcigar("MDDDM") == "1M3D1M"
    doAssert foldcigar("MMMMIDMMMM") == "4M1I1D4M"

  block:
    let sref = "CCATATGG"
    let squery = "CCAT**GG"
    let slen = cint(max(len(sref), len(squery)))
    let new_state_seq = newString(slen)# not newStringOfCap
    left_align_indels(sref, squery, slen, new_state_seq)
    doAssert new_state_seq == "MMDDMMMM"

  block:
    let sref = "CCAT**GG"
    let squery = "CCATATGG"
    let slen = cint(max(len(sref), len(squery)))
    let new_state_seq = newString(slen)# not newStringOfCap
    left_align_indels(sref, squery, 8, new_state_seq);
    doAssert new_state_seq == "MMIIMMMM"

  block:
    var sref = "CCATATGG*CC"
    let squery = "CCAT**GGGCC"
    let slen = cint(max(len(sref), len(squery)))
    let new_state_seq = newString(slen)# not newStringOfCap
    left_align_indels(sref, squery, slen, new_state_seq);
    doAssert new_state_seq == "MMDDMMIMMMM"

  block:
    # htsnim quals are offset already and seq[uint8].
    # the original viterbi expects char*.
    # found `cast [uint8_t](addr(bqs[0]))`
    # at https://forum.nim-lang.org/t/4647.
    # this is certainly faster then coding and decoding, but
    # required changes in c code, looks non-intuitive and depends on Nim's
    # internal seq representation.
    let def_qual: cint = 20
    let sref = "CCATATGG"
    let squery = "CCATGG"
    var bqual: seq[uint8]
    for i in 0..<len(squery):
      bqual.add(30)
    var alnseq = newString(max(len(sref), len(squery)))
    let shift = viterbi_c(sref, squery, cast[ptr uint8](addr(bqual[0])), alnseq, def_qual)
    doAssert shift == 0
    doAssert alnseq == "MMDDMMMM"

    echo "OK: all tests passed"