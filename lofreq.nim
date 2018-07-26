## LoFreq
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License


# read level filtering can only happen at the pileup stage.  so it
# makes sense to apply base level filtering there as well. in other
# words, here we only deal with filtered data. the data exchanged
# between plp and call is basically a quality histogram per event
# (base indel). the histogram allows a dense representation even for
# ultra high coverage data. the downside is, we cannot keep multiple
# qualities, because the linkage is broken in a histogram. so quality
# joining has to happen in the plp phase as well. as little filtering
# as possible should happen in the plp phase, because firstly LoFreq
# is meant to deal with noise and secondly you otherwise skew results
# (coverage etc.)

import tables
import cligen
import json
import parseutils
import math


##
## brief Computes log(exp(log_a) + exp(log_b))
##
## Taken from util.h of FAST source code:
## http://www.cs.cornell.edu/~keich/FAST/fast.tar.gz

proc log_sum(log_a: float32, log_b: float32): float32 =
  if log_a > log_b:
    log_a + ln(1.0 + exp(log_b-log_a))# FIXME use c log1p?
  else:
    log_b + ln(1.0 + exp(log_a-log_b))# FIXME use c log1p?


# FIXME: does it make sense to receive phred scores here instead?
# less conversion and no need to check range then?
proc pruned_prob_dist(err_probs: seq[float],# FIXME use ref to safe mem?
                      K: Natural,
                      bonf = 1.0, sig = 1.0): seq[float] =
  assert K > 0
  let N = err_probs.len
  var probvec = newSeq[float64](K+1)
  var probvec_prev = newSeq[float64](K+1)

  for f in err_probs:
    assert f>=0.0 and f<=1.0
  probvec_prev[0] = 0.0; # log(1.0)

  for n in 1..N:
    var k = 0;
    let pn = err_probs[n-1]
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
      probvec_prev[n] = -float(high(int))
      
    for k in countdown(min(n, K-1), 1):
      try:
        assert probvec_prev[k]<=0.0 and probvec_prev[k-1]<=0.0
      except AssertionError:
        raise
      probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                           probvec_prev[k-1] + log_pn)
      
    k = 0;
    assert probvec_prev[k] <= 0.0
    probvec[k] = probvec_prev[k] + log_1_pn
    
    if n==K:
      probvec[K] = probvec_prev[K-1] + log_pn;
      # FIXME prune here as well?
       
    elif n > K:
      try:
        assert probvec_prev[K] <= 0.0 and probvec_prev[K-1] <= 0.0 # FIXME
      except AssertionError:
        raise
      probvec[K] = log_sum(probvec_prev[K], probvec_prev[K-1] + log_pn)

      let pvalue = exp(probvec[K]);
          
      # FIXME store as phred scores instead?:
      # Q = -10*log_10(e^X), where X=probvec[K]
      # remember, log_b(x) = log_k(x)/log_k(b), i.e. log_10(Y) = log_e(Y)/log_e(10)
      # therefore, Q = -10 * log_e(e^X)/log_e(10) = -10 * X/log_e(10)
      # e.g.
      # >>> from math import log, log10, e
      # >>> X = -100
      # >>> -10 * log10(e**X)
      # 434.29448190325184
      # >>> -10 * X/log(10)
      # 434.2944819032518
      # */
      if pvalue * bonf > sig:
        return probvec# or nothing?
    swap(probvec, probvec_prev)

  # return prev because we just swapped (if not pruned)
  return probvec_prev


    
#Proc parse_plp_json(data: string) =
proc parse_plp_json(fname: string): Table[string, Table[int, int]]  =
  result = initTable[string, Table[int, int]]()
  
  # let dataJson = parseJson(data)
  let dataJson = parseFile(fname)
  
  # assert used in nim in action after parsing json string
  assert dataJson.kind == JObject

  for event, qhist in dataJson.pairs():
    for qual, count in qhist.pairs():
      var q = 0      
      assert qual[0] == 'Q'
      # FIXME: string slicing and toInt simply doesn't work. Try: var
      # x = readLine(stdin); let q = parseInt(x[1..len(x)])
      # parseSaturatedNatural works but looks incredibly awkward
      discard parseSaturatedNatural(qual, q, start=1)
      
      let c = count.getInt

      # FIXME allow for lowercase i.e bases mapping to rev strand.
      # Best done elsewhere and not here
      if not result.hasKey(event):
        result[event] = initTable[int, int]()
      assert result[event].hasKey(q) == false
      result[event][q] = c
    
  
proc main(plp_fname: string) =
  #var plpTable: Table[char, Table[int]]
  #for base in "ACGTN":
  #  plpTable[base] = initTable[int]()
    
  # populate with dummy data
  #for base in "ACGTN":
  #  for qual in [20, 30]:
  #    plpTable[base][qual] = int(qual/10)
  #echo(plpTable)
  # {'A': {20: 2, 30: 3}, 'C': {20: 2, 30: 3}, 'G': {20: 2, 30: 3}, 'N': {20: 2, 30: 3}, 'T': {20: 2, 30: 3}}
 

  #var data = %*
  #  [
  #     {"A": {"20": 2, "30": 3}, "C": {"20": 2, "30": 3}, "G": {"20": 2, "30": 3}, "N": {"20": 2, "30": 3}, "T": {"20": 2, "30": 3}, "+AC": {"20": 2, "30": 3}, "-TG": {"20": 2, "30": 3}}
  #  ]
  #parse_plp_json($data)
  # when passed as string we get a JArray?!
  # when passed as file we get JObject?!

  # plpTable =
  var plpTable = parse_plp_json(plp_fname)
  echo(plpTable)

  if false:
    var eprobs = @[0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    let num_failures = 1
    # produces 0.00995512 just like ppoibin
    let probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    let pvalue = exp(probvec[num_failures]);
    echo("num_failures=" & $num_failures & " pvalue=" & $pvalue)
  else:
    var eprobs = @[0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    for num_failures in countup(1, 10):
      let probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
      let pvalue = exp(probvec[num_failures]);
      echo("num_failures=" & $num_failures & " pvalue=" & $pvalue)
       
when isMainModule:
  import cligen
  dispatch(main)


