## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License

{.compile: "viterbi.c".}
{.passL: "-lm"}

# standard
import strformat
import tables

# third party
import hts

# project specific


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

    if rec.flag.secondary and skipSecondary:
      echo "DEBUG: secondary: write immediately or put in q"
      continue

    if rec.flag.unmapped:
      echo "DEBUG: unmapped: write immediately or put in q"
      continue

    var hasIndels = false
    for c in rec.cigar:
      if c.op in @[CigarOp.insert, CigarOp.deletion]:
        hasIndels = true
        break
    if not hasIndels:
      echo "DEBUG: noIndels write immediately or put in q"
      continue

    echo "DEBUG: shorten this mess by removing trailing soft clip and ignoring hard clips/padding"
    # make sure no internal soft, hard clip or padding
    # use slices to get bq and sq
    # keep start() and stop() positions as is
    # test on real data which has clips back
    # and move to func
    #type CigarOp* {.pure.} = enum match = 0'u32, insert, deletion, ref_skip, soft_clip, hard_clip, pad, equal, diff, back

    # remove soft clipped bases
    var queryOrig: string
    discard rec.sequence(queryOrig)
    var query = ""
    var bqualOrig: seq[uint8]
    discard rec.base_qualities(bqualOrig)
    var bqual: seq[uint8]

    var x = rec.start# coordinate on reference
    var y = 0;# coordinate on query

    # parse cigar string
    for c in rec.cigar:
      # or c.consumes.query and reference: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
      if c.op in @[CigarOp.match, CigarOp.equal, CigarOp.diff]:
        # FIXME translated from C. Can be shortened by slicing
        for j in 0..<c.len:
          query.add(queryOrig[y])
          bqual.add(bqualOrig[y])
          inc x
          inc y
      elif c.op == CigarOp.hard_clip:
        # in theory we should do nothing here but hard clipping info gets lost here FIXME
        echo "DEBUG: cannot handle hard clips yet"
        # FIXME write anyway FIXME next
      elif c.op == CigarOp.deletion:
        x += c.len
      elif c.op ==  CigarOp.insert:
        # FIXME translated from C. Can be shortened by slicing
        for j in 0..<c.len:
          query.add(queryOrig[y])
          bqual.add(bqualOrig[y])
          inc y
      elif c.op == CigarOp.soft_clip:
        y += c.len
      else:
        echo "DEBUG: unknown cigar op in read " & rec.qname
        # FIXME write anyway FIXME next

    echo "DEBUG: x=" & $x
    echo "DEBUG: y=" & $y
    echo "DEBUG: queryOrig=" & queryOrig
    echo "DEBUG: query=" & query
    echo "DEBUG: cigar=" & $rec.cigar

    if not refs.hasKey(chrom):
      echo "DEBUG: Loading " & chrom
      refs[chrom] = fai.get(chrom)

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
