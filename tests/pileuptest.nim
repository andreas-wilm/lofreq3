# standard library
import os
import json
import unittest
import strutils
import osproc
# project specific
import ../src/lofreqpkg/pileup/pileup
import ../src/lofreqpkg/pileup/processor
import ../src/lofreqpkg/pileup/postprocessing
# third party
import hts


suite "pileup tests":
  test "pileup functions":
    check mergeQuals(19, 20, 21) == 15

  test "test full pileup on samples":
    # this cannot be declared in the loop because it triggers a compiler bug
    var actual: JsonNode
    let lofreq = "../lofreq"

    for kind, folderPath in walkDir("pileup_samples"):
      actual = newJArray()
      echo "Testing " & folderPath
      let faiFile = joinPath(folderPath, "ref.fa")
      let bamFile = joinPath(folderPath, "alignments.bam")
      let expected = parseFile(joinPath(folderPath, "output.json"))

      # in the past we could call full_pileup with .then from pipelinetools
      # that stopped working with newer nim compilers, so we need to call
      # the binary. old: full_pileup(bamFile, faiFile, "ref:1-14", false)#.then(proc (x: JsonNode): void = actual.add(x)))

      let cmd = lofreq & " call -b " & bamfile & " -f " & faifile & " -r ref:1-13 --noMQ -p"
      var outp = execProcess(cmd, options = {poUsePath, poEvalCommand})# disables default poStdErrToStdOut
      outp = "[" & outp.replace("\n", ",") & "]"
      outp = outp.replace(",]", "]")
      let actual = parseJson(outp)
      check expected == actual
