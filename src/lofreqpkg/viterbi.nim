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


proc viterbi*(faFname: string, bamIn: string, bamOut: string) =
  var fai: Fai
  # keeping all observerd reference sequences in memory for speedup
  var refs = initTable[string, string]()
  var bam: Bam

  if not open(fai, faFname):
    quit("Could not open reference sequence file: " & faFname)

  open(bam, bamIn, fai=faFname)
  for rec in bam:
    var chrom = rec.chrom
    if not refs.hasKey(chrom):
      echo "DEBUG: Loading " & chrom
      refs[chrom] = fai.get(chrom)




#from lofreq_viterbi.c v2.1.4
#
#static int del_flag = 1;
#static int q2default = -1;
#static int reclip = 0;
#fetch_func(b, &tmp, del_flag, q2default, reclip);
#
#if read.unmapped:
#  write read
#  next
#
#if no indels in cigar:
#  write read
#  next
#
#if ref not loaded:
#  read reference and keep in mem for later
#extract reference with RWIN from read matchpos - RWIN to read endpos + RWIN
#
#FIXME what about hardclips?
#FIXME ignore secondary alignments?
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
