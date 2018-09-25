import unittest
import math

import ../src/lofreqpkg/call


# FIXME how to test main function in call (private and nameclash)?


suite "pvalue computation":

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
    
  setup:
    var pvalue: float
    var probvec: seq[float]
    var num_failures: int

    
  test "eprobs_10x30":
    var eprobs = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    
    num_failures = 1
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    check abs(pvalue - 0.00995512) < 1e-6

    num_failures = 2
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    check abs(pvalue - 4.476063e-05) < 1e-6

    
  test "eprobs_13-30":
    var eprobs = [0.050119, 0.039811, 0.031623, 0.025119, 0.019953, 0.015849, 0.012589, 0.010000, 0.007943, 0.006310, 0.005012, 0.003981, 0.003162, 0.002512, 0.001995, 0.001585, 0.001259, 0.001000,]

    num_failures = 1
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=1)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    check abs(pvalue - 0.2159726) < 1e-6

    num_failures = 2
    probvec = pruned_prob_dist(eprobs, num_failures, bonf=1.0, sig=0.05)
    pvalue = exp(probvec[num_failures]);
    #echo("DEBUG num_failures=" & $num_failures & " pvalue=" & $pvalue  & " prob2qual=" & $prob2qual(pvalue))
    check abs(pvalue - 0.02240387) < 1e-6
    
    
