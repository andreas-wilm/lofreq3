# standard library
import unittest
# third party
# /
# project specific
import  ../src/lofreqpkg/viterbi

suite "viterbi tests":
  test "Testing left-alignment of indels":
    var sref = "CCATATGG"
    var squery = "CCAT**GG"
    var slen = cint(max(len(sref), len(squery)))
    var new_state_seq = newString(slen)# not newStringOfCap
    left_align_indels(sref, squery, slen, new_state_seq)
    check new_state_seq == "MMDDMMMM"

    sref = "CCAT**GG"
    squery = "CCATATGG"
    slen = cint(max(len(sref), len(squery)))
    new_state_seq = newString(slen)# not newStringOfCap
    left_align_indels(sref, squery, 8, new_state_seq);
    check new_state_seq == "MMIIMMMM"

    sref =  "CCATATGG*CC"
    squery = "CCAT**GGGCC"
    slen = cint(max(len(sref), len(squery)))
    new_state_seq = newString(slen)# not newStringOfCap
    left_align_indels(sref, squery, slen, new_state_seq);
    check new_state_seq == "MMDDMMIMMMM"

#block:
#  echo "Testing viterbi realignment"
#  var def_qual: cint = 20
#  var sref = "CCATATGG"
#  var squery = "CCATGG"
#  var bqual = "??????"# = 30
#  var alnseq = newString(max(len(sref), len(squery)))
#  var shift = viterbi(sref, squery, bqual, alnseq, def_qual)
#  assert alnseq == "MMDDMMMM"

  test "Testing viterbi realignment":
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
    var shift = viterbi_c(sref, squery, cast[ptr uint8](addr(bqual[0])), alnseq, def_qual)
    check alnseq == "MMDDMMMM"

