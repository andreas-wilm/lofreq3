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
import processor
import recordFilter
import algorithm
import pipetools
import postprocessing
import times
import logging


#var L = newConsoleLogger(fmtStr = verboseFmtStr)
#addHandler(L)
# is there not way to set the fmtStr of the default handler?


proc full_pileup*(bamFname: string, faFname: string, ignBQ2: bool,
                  handler: DataToVoid) : void =
  ## Performs the pileup over all chromosomes listed in the bam file.
  var bam: Bam
  var fai: Fai

  if not open(fai, faFname):
    quit("Could not open fasta file.")

  if not open(bam, bamFname, index=true):
    quit("Could not open BAM file.")

  for chromosome in targets(bam.hdr):
    let name = chromosome.name
    var time: float

    var records = newRecordFilter(bam, name)
    
    time = cpuTime()
    var reference = fai.loadSequence(name)
    info("Time taken to load reference ", name, " ", cpuTime() - time)

    var storage = newSlidingDeque(name, handler)
    var processor = newProcessor(storage, ignBQ2)
    
    time = cpuTime()
    algorithm.pileup(records, reference, processor)
    info("Time taken to pileup reference ", name, " ", cpuTime() - time)

proc pileup*(bamFname: string, faFname: string, ignBQ2: bool = false) =
  full_pileup(bamFname, faFname, ignBQ2, toJson.then(print))


when isMainModule:
  import cligen
  dispatch(pileup)
