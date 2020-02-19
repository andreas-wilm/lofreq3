import strformat

{.compile: "viterbi.c".}
{.passL: "-lm"}

proc viterbi_test() {.cdecl, importc: "viterbi_test".}
# void
proc left_align_indels(sref: cstring, squery: cstring, slen: cint, new_state_seq: cstring) {.cdecl, importc: "left_align_indels".}
# returns shift
#proc viterbi*(sref: cstring, squery: cstring, bqual: cstring, saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}
proc viterbi(sref: cstring, squery: cstring, bqual: ptr uint8, saln: cstring, def_qual: cint): int {.cdecl, importc: "viterbi".}

block:
  echo "Testing left-alignment of indels"
  var sref = "CCATATGG"
  var squery = "CCAT**GG"
  var slen = cint(max(len(sref), len(squery)))
  var new_state_seq = newString(slen)# not newStringOfCap
  left_align_indels(sref, squery, slen, new_state_seq)
  assert new_state_seq == "MMDDMMMM"

  sref = "CCAT**GG"
  squery = "CCATATGG"
  slen = cint(max(len(sref), len(squery)))
  new_state_seq = newString(slen)# not newStringOfCap
  left_align_indels(sref, squery, 8, new_state_seq);
  assert new_state_seq == "MMIIMMMM"

  sref =  "CCATATGG*CC"
  squery = "CCAT**GGGCC"
  slen = cint(max(len(sref), len(squery)))
  new_state_seq = newString(slen)# not newStringOfCap
  left_align_indels(sref, squery, slen, new_state_seq);
  assert new_state_seq == "MMDDMMIMMMM"

#block:
#  echo "Testing viterbi realignment"
#  var def_qual: cint = 20
#  var sref = "CCATATGG"
#  var squery = "CCATGG"
#  var bqual = "??????"# = 30
#  var alnseq = newString(max(len(sref), len(squery)))
#  var shift = viterbi(sref, squery, bqual, alnseq, def_qual)
#  assert alnseq == "MMDDMMMM"

block:
  echo "Testing viterbi realignment"
  # htsnim quals are offset already and seq[uint8].
  # the original viterbi expects char*.
  # found `cast [uint8_t](addr(bqs[0]))`
  # at https://forum.nim-lang.org/t/4647.
  # this is certainly faster then coding and decoding, but
  # required changes in c code, looks non-intuitive and depends on Nim's
  # internal seq representation.
  var def_qual: cint = 20
  var sref = "CCATATGG"
  var squery = "CCATGG"
  var bqual: seq[uint8]
  for i in 0..<len(squery):
    bqual.add(30)
  var alnseq = newString(max(len(sref), len(squery)))
  var shift = viterbi(sref, squery, cast[ptr uint8](addr(bqual[0])), alnseq, def_qual)
  assert alnseq == "MMDDMMMM"


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
