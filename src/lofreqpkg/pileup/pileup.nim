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


proc full_pileup*(bamFname: string, faFname: string, ignBQ2: bool,
                  handler: DataToVoid) : void =
  ## Performs the pileup over all chromosomes listed in the bam file.
  var bam: Bam
  var fai: Fai

  if not open(fai, faFname):
    quit("Could not open fasta file.")

  open(bam, bamFname, index=true)

  for chromosome in targets(bam.hdr):
    let name = chromosome.name

    var records = newRecordFilter(bam, name)
    var reference = fai.loadSequence(name)

    var storage = newSlidingDeque(name, handler)
    var processor = newProcessor(storage, ignBQ2)

    algorithm.pileup(records, reference, processor)


proc pileup*(bamFname: string, faFname: string, ignBQ2: bool = false) =
  full_pileup(bamFname, faFname, ignBQ2, toJson.then(print))
 

when isMainModule:
  import cligen
  dispatch(pileup)
