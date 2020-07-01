## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License

# standard
import tables
import json
import math
import strutils
import logging
import strformat
# third party
from hts/stats import fishers_exact_test
# project specific
import vcf
import utils
import pileup/storage/containers/operationData
import pileup/storage/containers/qualityHistogram
import pileup/storage/containers/positionData

type VarType = enum snp, ins, del

const
  DEFAULT_MIN_VARQUAL* = 20
  DEFAULT_MIN_AF* = 0.005

type CallParams* = object
  minVarQual*: Natural
  minAF*: float


var callParams*: CallParams
callParams = CallParams(minVarQual:DEFAULT_MIN_VAR_QUAL,
                        minAF:DEFAULT_MIN_AF)

var logger = newConsoleLogger(fmtStr = verboseFmtStr,
                              useStderr = true)


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
  result = probVec[tailStartIndex]
  for i in tailStartIndex+1..<len(probVec):
    result = logSum(result, probVec[i]);


# FIXME add doc and refer to paper
# FIXME: does it make sense to receive phred scores here instead?
# less conversion and no need to check range then?
# FIXME no need to make public except for tests
# FIXME add pseudocode from paper here
# previously pruned_calc_prob_dist()
proc prunedProbDist(errProbs: openArray[float],# FIXME use ref to safe mem?
                     K: Natural, bonf = 1.0, sig = 1.0): seq[float] =
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
## and populate opsData using its set function
proc parseOperationData(node: JsonNode, opsData: var OperationData[string]) =
  for event, qHist in node.pairs():
    for qual, count in qHist.pairs():
      let q = parseInt($qual)
      let c = count.getInt
      opsData.set(event, q, c)


## brief parse pileup object from json
proc parsePlpJson(jsonString: string): PositionData =
  let dataJson = parseJson(jsonString)
  # the following assert is also used in 'nim in action' after parsing json string
  assert dataJson.kind == JObject

  let refIndex = dataJson["POS"].getInt
  let refBase = dataJson["REF"].getStr[0]
  let chromosome = dataJson["CHROM"].getStr
  result = newPositionData(refIndex, refBase, chromosome)

  parseOperationData(dataJson["M"], result.matches)
  parseOperationData(dataJson["I"], result.insertions)
  parseOperationData(dataJson["D"], result.deletions)

  # FIXME what about extra keys?
  # FIXME what about missing keys? try here?


proc setVarInfo(af: float, coverage: int, refBase: char, altBase: string,
  baseCountsStranded: CountTable[string], vtype: VarType): InfoField =
  result.af = af
  result.dp = coverage
  var dp4: Dp4
  dp4.refForward = baseCountsStranded.getOrDefault($refBase.toUpperAscii)
  dp4.refReverse = baseCountsStranded.getOrDefault($refBase.toLowerAscii)
  dp4.altForward = baseCountsStranded.getOrDefault(altBase.toUpperAscii)
  dp4.altReverse = baseCountsStranded.getOrDefault(altBase.toLowerAscii)
  if vtype == ins or vtype == del:# indel
    dp4.refForward = baseCountsStranded.getOrDefault($REF_SYMBOL_AT_INDEL_FW)
    dp4.refReverse = baseCountsStranded.getOrDefault($REF_SYMBOL_AT_INDEL_RV)
  result.dp4 = dp4
  var f = fishers_exact_test(dp4.refForward, dp4.refReverse,
    dp4.altForward, dp4.altReverse)
  result.sb = prob2qual(f.two)
  result.vtype = $vtype


proc getCountsAndEProbs[T](opData: T, vartype: VarType):
  (seq[float], Natural, CountTable[string], CountTable[string]) =
  # fill array of error probabilities, set baseCounts and baseCountsStranded.
  # note that base's strand is indicated by its case.
  var eProbs: seq[float] = @[]# base error probabilites
  var baseCounts = initCountTable[string]()# base counts
  var baseCountsStranded = initCountTable[string]()# strand aware counts
  var coverage: Natural = 0

  for base, qhist in pairs(opData.histogram):
    var thisBaseCount = 0
    if vartype == snp:
      assert len(base) == 1
    for qual, count in qhist:
      assert count>=0
      thisBaseCount += count
      # snp: '*' are deletions, i.e. physical coverage (count) with q=-1 (ignore)
      if vartype == snp and $base == $DEFAULT_BLANK_SYMBOL:
        continue
      assert qual>=0
      let e = qual2prob(qual)
      for i in countup(1, count):
        eProbs.add(e)
    baseCountsStranded.inc($base, thisBaseCount)
    baseCounts.inc($base.toUpperAscii, thisBaseCount)
    coverage += thisBaseCount

  return (eProbs, coverage, baseCounts, baseCountsStranded)


## result is a sequence, because we might return multiple variants for this position
proc callAtPos*(plp: PositionData): seq[Variant] =
  var eprobs: seq[float]
  var coverage: Natural
  var baseCounts: CountTable[string]
  var baseCountsStranded: CountTable[string]

  for vartype in low(VarType)..high(VarType):
    # FIXME there got to be an easier way to do this
    if vartype == snp:
      plp.matches.clean()
      (eprobs, coverage, baseCounts, baseCountsStranded) = getCountsAndEProbs(plp.matches, snp)
    elif vartype == ins:
      plp.insertions.clean()
      (eprobs, coverage, baseCounts, baseCountsStranded) = getCountsAndEProbs(plp.insertions, ins)
    elif vartype == del:
      plp.deletions.clean()
      (eprobs, coverage, baseCounts, baseCountsStranded) = getCountsAndEProbs(plp.deletions, del)
    else:
      raise newException(ValueError, "Illegal vartype" & $vartype)

    # determine valid alt bases and max alt count (not merged into above for readability)
    var altBases: seq[string]
    var maxAltCount = 0
    for b, c in pairs(baseCounts):
      if vartype == snp and (b[0] == plp.refBase or b[0] == 'N'):
        continue
      if b[0] in @[DEFAULT_BLANK_SYMBOL, REF_SYMBOL_AT_INDEL_FW, REF_SYMBOL_AT_INDEL_RV]:
        continue
      altBases.add(b)
      if c > maxAltCount:
        maxAltCount = c

    # loop over altBases and determine whether they are variants
    let maxAF = maxAltCount/coverage
    if maxAF >= callParams.minAF and maxAltCount > 0:# don't even compute probDist if we can't reach minAF with most abundant base
      logger.log(lvlDebug, fmt"Testing {vartype} at {plp.chromosome}:{plp.refIndex}: {baseCounts}")
      #logger.log(lvlDebug, fmt"eprobs  {eprobs}")
      let probVec = prunedProbDist(eProbs, maxAltCount)# FIXME call "by reference" to safe memory?
      var prevAltCount = high(int)# paranoid check to ensure sorting of pairs and early exit
      sort(baseCounts)
      for altBase, altCount in pairs(baseCounts):
        # weird: looping over altBases and getting altCount directly doesn't seem to work after
        # (destructive) sort on baseCounts. only iterators supported?
        if not altBases.contains(altBase):
          continue
        assert altCount <= prevAltCount
        prevAltCount = altCount

        # need to check minAF again, because test above was for most frequent base only
        let af = altCount/coverage
        if af < callParams.minAF or altCount == 0:
          #echo "DEBUG af<minAF for " & altBase & ":" & $altCount & " = " & $af & "<" & $minAF
          break# early exit possible since baseCounts are sorted

        # for maxAltCount exp(probVec[altCount]) == exp(probvecTailSum(probVec, altCount))
        let pvalue = exp(probvecTailSum(probVec, altCount))
        let qual = prob2qual(pvalue)
        logger.log(lvlDebug, fmt"af={af:.6f} altCount={altCount} for {vartype} {altBase} gives qual={qual}")
        if qual < callParams.minVarQual:
          #echo "DEBUG qual<minQual for " & altBase & ":" & $altCount & " = " & $qual & "<" & $minQual
          break# early exit possible since baseCounts are sorted

        var varRefBase, varAltBase: string
        if vartype == snp:
          varRefBase = $plp.refBase
          varAltBase = altBase
        elif vartype == ins:
          varRefBase = $plp.refBase
          varAltBase = $plp.refBase & altBase
        elif vartype == del:
          varRefBase = $plp.refBase & altBase
          varAltBase = $plp.refBase
        else:
          raise newException(ValueError, "Illegal vartype" & $vartype)

        var vcfVar = Variant(chrom : plp.chromosome, pos : plp.refIndex,
          id : ".", refBase : varRefBase, alt : varAltBase, qual : qual,
          filter : ".")
        vcfVar.info = setVarInfo(af, coverage, plp.refBase, altBase,
          baseCountsStranded, vartype)
        result.add(vcfVar)


proc call_from_plp*(plpFname: string, minVarQual: int = DEFAULT_MIN_VAR_QUAL,
                  minAF: float = DEFAULT_MIN_AF, logLevel = 0) =
  if logLevel >= 3:
    setLogFilter(lvlDebug)
  elif logLevel == 2:
    setLogFilter(lvlInfo)
  elif logLevel == 1:
    setLogFilter(lvlNotice)
  elif logLevel == 0:
    setLogFilter(lvlWarn)
  else:
    quit("Invalid log level")

  echo vcfHeader()

  callParams.minVarQual = minVarQual
  callParams.minAF = minAF

  var plpFh: File = if plpFname == "-": stdin else: open(plpFname)
  defer:
    if plpFh != stdin:
      plpFh.close

  for line in plpFh.lines:
    var plp = parsePlpJson(line)
    for v in callAtPos(plp):
      echo $v
  logger.log(lvlDebug, "Done. Goodbye")


when isMainModule:

  # pvalues were computed with [poibin](https://cran.r-project.org/package=poibin), e.g.
  # ```
  # > library("poibin")
  # > p = read.table('eprobs_10x30.txt')
  # > pp = c(p)$V1
  # nerrs = 1
  # pv = ppoibin(kk=length(pp)-nerrs, pp=1-pp)
  # pv
  # [1] 0.00995512
  # ```
  var pvalue: float
  var probvec: seq[float]
  var num_failures: int

  testblock "eprobs_10x30":
    var eprobs = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

    num_failures = 1
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    doAssert abs(pvalue - 0.00995512) < 1e-6

    num_failures = 2
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    doAssert abs(pvalue - 4.476063e-05) < 1e-6

  testblock "eprobs_13-30":
    var eprobs = [0.050119, 0.039811, 0.031623, 0.025119, 0.019953, 0.015849, 0.012589, 0.010000, 0.007943, 0.006310, 0.005012, 0.003981, 0.003162, 0.002512, 0.001995, 0.001585, 0.001259, 0.001000,]

    num_failures = 1
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=1)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    doAssert abs(pvalue - 0.2159726) < 1e-6

    num_failures = 2
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    doAssert abs(pvalue - 0.02240387) < 1e-6
  
  echo "OK: all tests passed"