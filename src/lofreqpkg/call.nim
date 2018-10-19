## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License


import tables
import json
import utils
import math
import strutils

import vcf


type QualHist = Table[string, CountTable[int]]


type PositionData = object
    ## The 'PositionData' object keeping all information concerning one parti-
    ## cular position on the reference.
    refIndex: int
    refBase: char
    chromosome: string
    matches: QualHist
    deletions: QualHist
    insertions: QualHist


## brief Computes log(exp(logA) + exp(logB))
##
## Taken from util.h of FAST source code:
## http://www.cs.cornell.edu/~keich/FAST/fast.tar.gz
proc logSum(logA: float32, logB: float32): float32 =
  if logA > logB:
    logA + ln(1.0 + exp(logB-logA))# FIXME use c log1p?
  else:
    logB + ln(1.0 + exp(logA-logB))# FIXME use c log1p?


# FIXME add doc and refer to paper
# FIXME: does it make sense to receive phred scores here instead?
# less conversion and no need to check range then?
# FIXME no need to make public except for tests
# FIXME add pseudocode from paper here
# previously pruned_calc_prob_dist()
proc prunedProbDist*(errProbs: openArray[float],# FIXME use ref to safe mem?
                     K: Natural,
                     bonf = 1.0, sig = 1.0): seq[float] =
  assert K > 0
  var probVec = newSeq[float64](K+1)
  var probVecPrev = newSeq[float64](K+1)

  for f in errProbs:
    assert f>=0.0 and f<=1.0
  probVecPrev[0] = 0.0; # log(1.0)

  for n in 1..errProbs.len:
    var k = 0;
    let pn = errProbs[n-1]
    var log_pn = 0.0'f64
    var log_1_pn = 0.0'f64

    # if pn=0, log() will fail (we shouldn't see pn=0 though, because
    # phred scores). likewise if pn=1 (Q0!) then log1p(-pn) = log(1-1)
    # = log(0) will fail. therefore, test
    if pn == 0.0:
       log_pn = -float(high(int))# log(DBL_EPSILON)
    else:
       log_pn = ln(pn)
    if pn == 1.0:
       log_1_pn = -float(high(int))# log_1_pn = log1p(-pn+DBL_EPSILON);
    else:
       log_1_pn = ln(1.0-pn)# /* 0.0 = log(1.0) # FIXME use c log1p?

    if n < K:
      probVecPrev[n] = -float(high(int))

    for k in countdown(min(n, K-1), 1):
      assert probVecPrev[k]<=0.0 and probVecPrev[k-1]<=0.0
      probvec[k] = logSum(probVecPrev[k] + log_1_pn,
                           probVecPrev[k-1] + log_pn)

    k = 0;
    assert probVecPrev[k] <= 0.0
    probvec[k] = probVecPrev[k] + log_1_pn

    if n==K:
      probvec[K] = probVecPrev[K-1] + log_pn;
      # FIXME prune here as well?

    elif n > K:
      assert probVecPrev[K] <= 0.0 and probVecPrev[K-1] <= 0.0
      probvec[K] = logSum(probVecPrev[K], probVecPrev[K-1] + log_pn)

      let pvalue = exp(probvec[K]);

      # FIXME store as phred scores instead?
      # Q = -10*log_10(e^X), where X=probvec[K]
      # remember, logB(x) = log_k(x)/log_k(b), i.e. log_10(Y) = log_e(Y)/log_e(10)
      # therefore, Q = -10 * log_e(e^X)/log_e(10) = -10 * X/log_e(10)
      # e.g.
      # >>> from math import log, log10, e
      # >>> X = -100
      # >>> -10 * log10(e**X)
      # 434.29448190325184
      # >>> -10 * X/log(10)
      # 434.2944819032518

      # early exit
      if pvalue * bonf > sig:
        # explicitly limiting to valid range
        return probvec[0..K]
    swap(probvec, probVecPrev)

  # return prev because we just swapped (if not pruned)
  return probVecPrev[0..K] # explicitly limiting to valid range


proc parseOperationData(node: JsonNode): QualHist =
  # FIXME parse reverse and forward counts here
  # FIXME wrongly encoded in json
  result = initTable[string, CountTable[int]]()

  for event, qHist in node.pairs():
    discard result.hasKeyOrPut(event, initCountTable[int]())
    for qual, count in qHist.pairs():
      ##assert qual[0] == 'Q'
      ## FIXME: string slicing and toInt simply doesn't work. Try: var
      ## x = readLine(stdin); let q = parseInt(x[1..len(x)])
      ## parseSaturatedNatural works but looks incredibly awkward
      #discard parseSaturatedNatural(qual, q, start=1)
      let q = parseInt($qual)
      let c = count.getInt

      assert result[event].hasKey(q) == false
      result[event][q] = c
  #echo("DEBUG returning " & $result)


## brief parse pileup object from json
proc parsePlpJson(jsonString: string): PositionData =
  let dataJson = parseJson(jsonString)
  #let dataJson = parseFile(fname)

  # assert used in nim in action after parsing json string
  assert dataJson.kind == JObject
  # FIXME should we validate the keys before starting to parse?
  # what about extra keys for example. should they be ignored?

  result.chromosome = dataJson["chromosome"].getStr
  result.refIndex =  dataJson["referenceIndex"].getInt
  result.refBase = dataJson["referenceBase"].getStr[0].toUpperAscii()

  result.matches = parseOperationData(dataJson["matches"])
  result.insertions = parseOperationData(dataJson["insertions"])
  result.deletions = parseOperationData(dataJson["deletions"])


proc call*(plpFname: string, minQual: int = 20, minAF: float = 0.005) =
  echo vcfHeader()

  for line in lines(plpFname):
    var plp = parsePlpJson(line)
    var baseCounts = initTable[string, int]()
    var baseCount = 0
    var eProbs: seq[float] = @[]
    var coverage = 0

    # fill array of error probabilities and set baseCounts[] as well as coverage
    for base, qhist in plp.matches.pairs():
      assert len(base) == 1# because SNP branch
      baseCount = 0
      for qual, count in qhist:
        assert count>=0
        baseCount += count
        # '*' are deletions, i.e. physical coverage (count) with q=-1 (ignore)
        if base != "*":
          assert qual>=0
          let e = qual2prob(qual)
          for i in countup(1, count):
            eProbs.add(e)
      let uBase = base.toUpperAscii
      discard baseCounts.hasKeyOrPut(uBase, 0)
      baseCounts[uBase] += baseCount
      coverage += baseCount
    assert len(eProbs) + baseCounts.getOrDefault("*") == coverage

    # determine max number of alt snp counts. used for pruning in prunedProbDist()
    var maxAltCount = 0
    for b, c in pairs baseCounts:
      assert len(b) == 1
      if b[0] != plp.refBase and b[0] != 'N' and b[0] != '*':
        if c > maxAltCount:
          maxAltCount = c

    var pvalue = -1.0
    var varQual = -1
    if maxAltCount > 0:# chance to exit early if minAF is not met
      # first pass loop: skip probDist calculation if below minAF
      var passMinAF = false
      for b, c in pairs baseCounts:
        if b[0] != plp.refBase and b[0] != 'N' and b[0] != '*':
          let af = baseCounts[b]/coverage
          if af >= minAF:
            passMinAF = true
            break
      if passMinAF == true:
        # FIXME call by reference to safe memory
        let probVec = prunedProbDist(eProbs, maxAltCount)
        for b, c in pairs baseCounts:
          if b[0] != plp.refBase and b[0] != 'N' and b[0] != '*':
            # need to check minAF again, because checks above are for any base in this column
            let af = baseCounts[b]/coverage
            if af >= minAF:
              let pvalue = exp(probVec[maxAltCount]);
              let qual = prob2qual(pvalue)
              if qual >= minQual:
                var vcfVar: Variant
                vcfVar.chrom = plp.chromosome
                vcfVar.pos = plp.refIndex
                vcfVar.id = "."
                vcfVar.refBase = plp.refBase
                vcfVar.alt = b
                vcfVar.qual = qual
                vcfVar.filter = "."
                var info: InfoField
                info.af = af
                info.sb = 0# FIXME
                var dp4: Dp4
                dp4.refForward = 0#FIXME
                dp4.refReverse = 0#FIXME
                dp4.altForward = 0#FIXME
                dp4.altReverse = 0#FIXME
                info.dp4 = dp4
                info.indel = false
                vcfVar.info = info
                echo "DEBUG: ", plp.chromosome, " ", plp.refIndex, " ", plp.refBase, " ", coverage, " ", baseCounts, " ", maxAltCount, " ", len(eProbs), " ", pvalue, " ", varQual
                echo $vcfVar
  

when isMainModule:
  import cligen
  dispatch(call)
