## Implementation of an alternative procedure to samtools' pileup that does not
## require keeping pointers to all the reads in memory and provides support for
## streaming results to microservices as soon as possible. This was inspired by
## https://brentp.github.io/post/no-pile/.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

import os
import hts
import interfaces/iSequence
import storage/slidingDeque
import recordFilter
import algorithm
import pipetools
import postprocessing
import times
import logging


#var L = newConsoleLogger(fmtStr = verboseFmtStr)
#addHandler(L)
# is there not way to set the fmtStr of the default handler?


proc full_pileup*(bamFname: string, faFname: string, handler: DataToVoid) : void =
  ## Performs the pileup over all chromosomes listed in the bam file.
  var bam: Bam
  var fai: Fai

  if not open(bam, bamFname, index=true):
    quit("Could not open BAM file " & bamFname)

  if not open(fai, faFname):
    quit("Could not open reference " & faFname)

  # FIXME in-lieu of region
  for chromosome in targets(bam.hdr):
    var records = newRecordFilter(bam, chromosome.name)

    let time = cpuTime()
    algorithm.pileup(fai, records, handler)
    info("Time taken to pileup reference ", chromosome.name, " ", cpuTime() - time)


proc pileup*(bamFname: string, faFname: string) =
  #full_pileup(bamFname, faFname, doNothing)
  full_pileup(bamFname, faFname, toJson.then(print))


when isMainModule:
  import cligen
  dispatch(pileup)
