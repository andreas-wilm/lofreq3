## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License


import tables
import json
import parseutils
import utils
import math


## brief pileup object parsed from json
type PlpObj = object
  sq: string
  pos: Natural
  refStr: string
  qHist: Table[string, Table[int, int]]


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
      try:# FIXME why try?
        assert probVecPrev[k]<=0.0 and probVecPrev[k-1]<=0.0
      except AssertionError:
        raise
      probvec[k] = logSum(probVecPrev[k] + log_1_pn,
                           probVecPrev[k-1] + log_pn)
      
    k = 0;
    assert probVecPrev[k] <= 0.0
    probvec[k] = probVecPrev[k] + log_1_pn
    
    if n==K:
      probvec[K] = probVecPrev[K-1] + log_pn;
      # FIXME prune here as well?
       
    elif n > K:
      try:# FIXME why try?
        assert probVecPrev[K] <= 0.0 and probVecPrev[K-1] <= 0.0 # FIXME
      except AssertionError:
        raise
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


## brief parse pileup object from json
proc parsePlpJson(fname: string): PlpObj =
  result.qHist = initTable[string, Table[int, int]]()
  
  # let dataJson = parseJson(data)
  let dataJson = parseFile(fname)
  
  # assert used in nim in action after parsing json string
  assert dataJson.kind == JObject
  # FIXME should we validate the keys before starting to parse?
  # what about extra keys for example. should they be ignored?
  
  result.sq = dataJson["sq"].getStr
  result.pos =  dataJson["pos"].getInt
  result.refStr = dataJson["ref"].getStr
  
  for event, qHist in dataJson["qhist"].pairs():
    for qual, count in qHist.pairs():
      var q = 0      
      assert qual[0] == 'Q'
      # FIXME: string slicing and toInt simply doesn't work. Try: var
      # x = readLine(stdin); let q = parseInt(x[1..len(x)])
      # parseSaturatedNatural works but looks incredibly awkward
      discard parseSaturatedNatural(qual, q, start=1)
      
      let c = count.getInt

      # FIXME allow for lowercase i.e bases mapping to rev strand.
      # Best done elsewhere and not here
      if not result.qHist.hasKey(event):
        result.qHist[event] = initTable[int, int]()
      assert result.qHist[event].hasKey(q) == false
      result.qHist[event][q] = c


## brief create vcf header
# FIXME check conformity. refFa and src likely missing since only known in plp
proc vcfHeader(src: string = "", refFa: string = ""): string =
  result = "##fileformat=VCFv4.2\n"
  result = result & "##fileDate=" & dateStr() & "\n"
  if len(src)>0:
    result = result & "##source=" & src & "\n"
  if len(refFa)>0:
    result = result & "##reference=" & refFa & "\n"

  result = result & """##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO""" 

    
proc call*(plpFname: string) =
  echo vcfHeader()
  var plp = parsePlpJson(plpFname)
  echo("FIXME: compute pvalue for ", plp)
  
  
when isMainModule:
  import cligen
  dispatch(call)
