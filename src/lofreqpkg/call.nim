## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License

# standard
import tables
import json
import utils
import math
import strutils
# third party
from hts/stats import fishers_exact_test
# project specific
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


## brief Computes sum of probvec values (log space) starting from (including)
## tail_startindex to (excluding) probvec_len
##
proc probvecTailSum(probVec: openArray[float], tailStartIndex: int): float =
  result = probVec[tailStartIndex];
  for i in tailStartIndex+1..<len(probVec):
    result = logSum(result, probVec[i]);


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
  # explicitly limiting to valid range
  return probVecPrev[0..K]


## brief parse quality histogram for all operations from json node
proc parseOperationData(node: JsonNode): QualHist =
  result = initTable[string, CountTable[int]]()

  for event, qHist in node.pairs():
    discard result.hasKeyOrPut(event, initCountTable[int]())
    for qual, count in qHist.pairs():
      let q = parseInt($qual)
      let c = count.getInt

      assert result[event].hasKey(q) == false
      result[event][q] = c


## brief parse pileup object from json
proc parsePlpJson(jsonString: string): PositionData =
  let dataJson = parseJson(jsonString)
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
    var baseCounts = initCountTable[string]()# base counts
    var baseCountsStranded = initCountTable[string]()# strand aware counts
    var eProbs: seq[float] = @[]
    var coverage = 0

    # Fill array of error probabilities, set baseCounts, baseCountsStranded and coverage.
    # Note that base's strand is indicated by its case.
    for base, qhist in plp.matches.pairs():
      var thisBaseCount = 0
      assert len(base) == 1# because SNP branch
      for qual, count in qhist:
        assert count>=0
        thisBaseCount += count
        # '*' are deletions, i.e. physical coverage (count) with q=-1 (ignore)
        if base != "*":
          assert qual>=0# "*" have q=-1
          let e = qual2prob(qual)
          for i in countup(1, count):
            eProbs.add(e)
      coverage += thisBaseCount
      baseCountsStranded.inc(base, thisBaseCount)
      baseCounts.inc(base.toUpperAscii, thisBaseCount)
    assert len(eProbs) + baseCounts.getOrDefault("*") == coverage

    # determine max number of alt snp counts. used for pruning in prunedProbDist().
    let (maxAltBase, maxAltCount) = largest(baseCounts)
    if maxAltCount > 0:
      # skip probDist calculation if below minAF
      var passMinAF = false
      for b, c in pairs baseCounts:
        if b[0] != plp.refBase and b[0] != 'N' and b[0] != '*':# FIXME code duplication
          let af = baseCounts[b]/coverage
          if af >= minAF:
            passMinAF = true
            break
      if passMinAF == true:
        # FIXME call "by reference" to safe memory?
        let probVec = prunedProbDist(eProbs, maxAltCount)
        var prevAltCount = high(int)# paranoid check to ensure sorting of pairs and early exit
        sort(baseCounts)
        for altBase, altCount in pairs baseCounts:
          if altBase[0] != plp.refBase and altBase[0] != 'N' and altBase[0] != '*':# FIXME code duplication
            assert altCount <= prevAltCount
            prevAltCount = altCount
            # need to check minAF again, because test above was for more frequent base only
            let af = altCount/coverage
            if af >= minAF:
              #let pvalue = exp(probVec[altCount]);
              let pvalue = exp(probvecTailSum(probVec, altCount))
              let qual = prob2qual(pvalue)
              if qual >= minQual:
                var vcfVar = Variant(chrom : plp.chromosome, pos : plp.refIndex,
                  id : ".", refBase : plp.refBase, alt : altBase, qual : qual, filter : ".")
                var info: InfoField
                info.af = af
                info.dp = coverage
                var dp4: Dp4
                dp4.refForward = baseCountsStranded.getOrDefault($plp.refBase.toUpperAscii)
                dp4.refReverse = baseCountsStranded.getOrDefault($plp.refBase.toLowerAscii)
                dp4.altForward = baseCountsStranded.getOrDefault(altBase.toUpperAscii)
                dp4.altReverse = baseCountsStranded.getOrDefault(altBase.toLowerAscii)
                info.dp4 = dp4

                info.sb = 0# FIXME against all or only DP4?
                var f = fishers_exact_test(dp4.refForward, dp4.refReverse, dp4.altForward, dp4.altReverse)
                info.sb = prob2qual(f.two)
                vcfVar.info = info
                echo $vcfVar
              else:
                # early exit possible since baseCounts sorted
                break
            else:
              # early exit possible since baseCounts sorted
              break

when isMainModule:
  import cligen
  dispatch(call)
