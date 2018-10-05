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


proc pileup*(bamFname: string, faFname: string) =
  ## Performs the pileup over all chromosomes listed in the bam file.
  ## FIXME: enable multithreading.
  var bam: Bam
  var fai: Fai

  if not open(fai, faFname):
    quit("Could not open fasta file.")
  
  open(bam, bamFname, index=true)

  for chromosome in targets(bam.hdr):
    let name = chromosome.name

    var records = newRecordFilter(bam, name)
    var reference = fai.loadSequence(name)
    
    var storage = newSlidingDeque(name, toJson.then(print))
    var processor = newProcessor(storage)
    pileup(records, reference, processor)

when isMainModule:
  import cligen
  dispatch(pileup)
