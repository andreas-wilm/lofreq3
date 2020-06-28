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

# third party
import hts

# project specific
#/

# void
proc left_align_indels*(sref: cstring, squery: cstring, slen: cint, new_state_seq: cstring) {.cdecl, importc: "left_align_indels".}


# returns shift
#proc viterbi*(sref: cstring, squery: cstring, bqual: cstring, saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}
proc viterbi_c*(sref: cstring, squery: cstring, bqual: ptr uint8, saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}


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
# FIXME add tests
proc getRefContext(rec: Record, refSq: string, refPadding: int): string =
  var s = rec.start - refPadding
  if s < 0:
    s = 0
  var e = rec.stop + refPadding
  if e >= len(refSq):
    e = len(refSq)
  result = toUpperAscii(refSq[s..e])#toupper to avoid ref masking


proc encCigarChar(op: char): uint32 =
  let idx = "MIDNSHP=XB".find(op)
  doAssert idx >= 0
  result = uint32(idx)


proc foldCigar(rawCigar: string): seq[CigarElement] =
  # translation of LoFreq2 function
  var endIdx = len(rawCigar)-1# this is only needed because the rawCigar is overallocated and len doesn't work!?
  for i, opChar in rawCigar.pairs:
    if not opChar.isUpperAscii:
      endIdx = i-1
      break
    doAssert opChar in "MID"# viterbi should not have produced anything else

  var curr_op = encCigarChar(rawCigar[0])
  var curr_oplen = 1
  for op_char in rawCigar[1..endIdx]:
    let this_op = encCigarChar(rawCigar[0])
    if this_op != curr_op:
      result.add(CigarElement(uint32(curr_oplen shl 4) or curr_op))
      curr_op = this_op
      curr_oplen = 1
    else:
      curr_oplen += 1
  result.add(CigarElement(uint32(curr_oplen shl 4) or curr_op))


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
    #echo "DEBUG cigar=" & $rec.cigar & " leadingSkipOps=" & $leadingSkipOps & " trailingSkipOps=" & $trailingSkipOps

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
    #echo "DEBUG leadingSoftClipLen=" & $leadingSoftClipLen & " trailingSoftClipLen=" & $trailingSoftClipLen

    return (leadingSkipOps, leadingSoftClipLen, 
      trailingSkipOps, trailingSoftClipLen)


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
      #echo "DEBUG: Loading " & chrom
      refs[chrom] = fai.get(chrom)

    stderr.writeLine("DEBUG incoming read = " & $rec.tostring())

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

    let q2def: cint = 2
    stderr.writeLine("WARN: fixed q2def=" & $q2def)
    stderr.writeLine("WARN: implement median for testing.")
    let realnCigarRaw = newString(max(len(queryWOSoftClip), len(refContext)))
    # WARNING with newStringOfCap we get a segfault, but now newCigarRaw is oversized
    # printing works, but len is incorrect. creating a copy to a second string with $ doesn't help
    # FIXME what to do with shift?
    let shift = viterbi_c(refContext, queryWOSoftClip,
      cast[ptr uint8](addr(bqualWOSoftCLip[0])), realnCigarRaw, q2def)
    #echo "DEBUG realnCigarRaw ='" & $realnCigarRaw & "'"

    let realnCigar = foldCigar(realnCigarRaw)
    #echo "DEBUG realnCigar=" & $realnCigarRaw
    var fullRealnCigarSeq = concat(leadingSkipOps, realnCigar, trailingSkipOps)
    # https://github.com/brentp/hts-nim/commit/0bf2683b17f5ccc27a7af4731f187452b287cb61
    GC_ref(fullRealnCigarSeq)
    let fullRealnCigar = newCigar(fullRealnCigarSeq)
    #echo "DEBUG fullRealnCigar=" & $fullRealnCigar

    # check if start pos was shifted
    var realnStart = rec.start
    var lower = rec.start - refPadding;
    if lower < 0:
      lower = 0
    if shift - (rec.start - lower) != 0:
      realnStart = rec.start + (shift - (rec.start - lower))
      stderr.writeLine("DEBUG: new realnStart = " & $realnStart & " lower=" & $lower & " shift=" & $shift)
    else:
      stderr.writeLine("DEBUG: new=old realnStart = " & $realnStart)
      
    echo createRealnRec(rec, realnStart, fullRealnCigar)

    #oBam.close()

    GC_unref(fullRealnCigarSeq)


stderr.writeLine("FIXME test case: shift of read")
stderr.writeLine("FIXME test case: skip and skipS at both sides")
stderr.writeLine("WARN wrongly changed offset in ./lofreq viterbi -f tests/viterbi/NC_011770_head.fa -i tests/viterbi/pseudomonas_pair_screwed_up_cigar.bam")


when isMainModule:
  import cligen
  dispatch(viterbi)
