import os
import json
import hts
import unittest

import ../src/lofreqpkg/pileup/pileup
import ../src/lofreqpkg/pileup/postprocessing
import ../src/lofreqpkg/pileup/pipetools



suite "pileup tests":

  test "test full pileup on samples":
    # this cannot be declared in the loop because it triggers a compiler bug
    var actual: JsonNode

    for kind, folderPath in walkDir("pileup_samples"):
      actual = newJArray()

      let faiFile = joinPath(folderPath, "ref.fa")
      let bamFile = joinPath(folderPath, "alignments.bam")

      let expected = parseFile(joinPath(folderPath, "output.json"))

      full_pileup(bamFile, faiFile, false,
                  toJson.then(proc (x: JsonNode): void = actual.add(x)))

      check expected == actual

