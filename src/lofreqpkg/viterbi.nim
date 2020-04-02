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

# third party
import hts

# project specific
#/

# void
proc left_align_indels*(sref: cstring, squery: cstring, slen: cint, new_state_seq: cstring) {.cdecl, importc: "left_align_indels".}


# returns shift
#proc viterbi*(sref: cstring, squery: cstring, bqual: cstring, saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}
proc viterbi_c*(sref: cstring, squery: cstring, bqual: ptr uint8, saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}


proc viterbi*(faFname: string, bamIn: string, bamOut: string, skipSecondary = true) =
  var fai: Fai
  # keeping all observerd reference sequences in memory for speedup
  var refs = initTable[string, string]()
  var bam: Bam

  if not open(fai, faFname):
    quit("Could not open reference sequence file: " & faFname)

  open(bam, bamIn, fai=faFname)
  for rec in bam:
    var chrom = rec.chrom

    # skip secondary reads unless requested
    if rec.flag.secondary and skipSecondary:
      echo "DEBUG: secondary: write immediately or put in q"
      continue

    # skip unmapped reads
    if rec.flag.unmapped:
      echo "DEBUG: unmapped: write immediately or put in q"
      continue

    # skip reads without indels
    var hasIndels = false
    for c in rec.cigar:
      if c.op in @[CigarOp.insert, CigarOp.deletion]:
        hasIndels = true
        break
    if not hasIndels:
      echo "DEBUG: noIndels write immediately or put in q"
      continue

    # load reference if not cached
    if not refs.hasKey(chrom):
      echo "DEBUG: Loading " & chrom
      refs[chrom] = fai.get(chrom)

    # get reference context for aligned read
    let refPadding = 10
    var s = rec.start - refPadding
    if s < 0:
      s = 0
    var e = rec.stop + refPadding
    if e >= len(refs[chrom]):
      e = len(refs[chrom])
    var refContext = toUpperAscii(refs[chrom][s..e])#toupper to avoid ref masking
    echo "refContext=" & refContext

    # make sure no internal soft, hard clip or padding
    # use slices to get bq and sq
    # keep start() and stop() positions as is
    # test on real data which has clips back
    # and move to func
    #type CigarOp* {.pure.} = enum match = 0'u32, insert, deletion, ref_skip, soft_clip, hard_clip, pad, equal, diff, back

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
    echo "DEBUG cigar=" & $rec.cigar & " leadingSkipOps=" & $leadingSkipOps & " trailingSkipOps=" & $trailingSkipOps

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
    echo "DEBUG leadingSoftClipLen=" & $leadingSoftClipLen & " trailingSoftClipLen=" & $trailingSoftClipLen

    var query: string
    discard rec.sequence(query)
    var queryWOSoftClip = ""
    queryWOSoftClip = query[leadingSoftClipLen .. len(query)-1-trailingSoftClipLen]

    var bqual: seq[uint8]
    discard rec.base_qualities(bqual)
    var bqualWOSoftCLip: seq[uint8]
    bqualWOSoftCLip = bqual[leadingSoftClipLen .. len(bqual)-1-trailingSoftClipLen]

    echo "DEBUG: origQuery=" & query
    echo "DEBUG: queryWOSoftClip=" & queryWOSoftClip# & " bqualWOSoftCLip=" & $bqualWOSoftCLip

    let q2def: cint = 2
    echo "WARN: fixed q2def=" & $q2def
    echo "WARN: implement median for testing."
    var newCigarRaw = newString(max(len(queryWOSoftClip), len(refContext)))

    var shift = viterbi_c(refContext, queryWOSoftClip, cast[ptr uint8](addr(bqualWOSoftCLip[0])), newCigarRaw, q2def)
    echo "DEBUG newCigarRaw=" & $newCigarRaw


    echo "NotImplemented"
    # reconstruct cigar: take care of padding, hardclip and softclip
  # if new cigar != old cigar
  echo "WARNING: mapper generated scores will be invalid / off (NM, MC, MD, AS and even MQ) for realigned reads"
# ../lofreq.git/src/lofreq/lofreq_viterbi.c
#from lofreq_viterbi.c v2.1.4
#
#determine read-refstartpos ignoring clips
#determine read-refendpos ignoring clips
#
#q2def = user arg or median of all bases
#int shift = viterbi(ref, query, bqual, aln, q2def);
#if shift:
#  update start read-refstartpos
#compute cigar:
#  add softclip at start
#  add viterbi-cigar
#  add softclip at end
#update cigar
#
#if newcigar != oldcigar:
#  for tag in NM MC MD AS
#    delete tag
#    # FIXME are there more?
#
#write updated read


when isMainModule:
  import cligen
  dispatch(viterbi)
