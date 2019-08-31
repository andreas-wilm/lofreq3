## Implementation of an alternative procedure to samtools' pileup that does not
## require keeping pointers to all the reads in memory and provides support for
## streaming results to microservices as soon as possible. This was inspired by
## https://brentp.github.io/post/no-pile/.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

# standard
import os
import hts
import interfaces/iSequence
import storage/slidingDeque
import recordFilter
import algorithm
import postprocessing
import times
import logging
import strutils

# project specific
import ../region


var logger = newConsoleLogger(fmtStr = verboseFmtStr, useStderr = true)


proc auto_fill_region*(reg: Region, targets: seq[Target]): Region = 
  result = reg
  var foundTarget = false
  for t in targets:
    if t.name == reg.sq:
      result.s = 0
      result.e = (int)t.length
      foundTarget = true
      break
  if not foundTarget:
    raise newException(ValueError, "Couldn't find " & reg.sq & " in BAM header. Valid entries are " & $targets)


proc full_pileup*(bamFname: string, regions = "", faFname = "", useMQ: bool, handler: DataToVoid) : void =
  ## Performs the pileup over all chromosomes listed in the bam file.
  var bam: Bam
  var fai: Fai

  if not open(bam, bamFname, index=true):
    quit("Could not open BAM file " & bamFname)

  if len(faFname)!=0:
    logger.log(lvlInfo, "Opening index for " & faFname)
    if not open(fai, faFname):
      quit("Could not open reference " & faFname)
  else:
    logger.log(lvlInfo, "No reference file given")

  for regstr in regions.split(','):
    var reg = reg_from_str(reg_str)

    if reg.s == 0 and reg.e == 0:# only sq given instead of full region
      var targets = targets(bam.hdr)
      reg = auto_fill_region(reg, targets)
    logger.log(lvlInfo, "Starting pileup for " & $reg)

    var records = newRecordFilter(bam, reg.sq, reg.s, reg.e)

    let time = cpuTime()
    algorithm.pileup(fai, records, reg, useMQ, handler)
    info("Time taken to pileup reference ", reg.sq, " ", cpuTime() - time)


proc pileup*(bamFname: string, regions: string, noMQ: bool = false, faFname = "") =
  #full_pileup(bamFname, faFname, doNothing)
  full_pileup(bamFname, regions, faFname, not noMQ, toJsonAndPrint)


when isMainModule:
  import cligen
  dispatch(pileup)
