## Implementation of an alternative procedure to samtools' pileup that does not require 
## keeping pointers to all the reads in memory and provides support for streaming results
## to microservices as soon as possible. This was inspired by https://brentp.github.io/post/no-pile/.
import os
import messaging

import hts
import interfaces/iSequence
import storage/containers/positionData
import storage/slidingDeque
import recordFilter
import jsonCollector
import algorithm

proc pileup*(bamFname: string, faFname: string) =
  var bam: Bam
  var fai: Fai

  if not open(fai, faFname):
    quit("Could not open fasta file.")
  
  open(bam, bamFname, index=true)
  
  for chromosome in targets(bam.hdr):
    let name = chromosome.name
  
    var records = newRecordFilter(bam, name)
    var reference = fai.getISequence(name)
    var storage = newSlidingDeque(200, newJsonCollector(name).getICollector)
  
    pileup(records, reference, storage)

    
when isMainModule:
  import cligen
  dispatch(pileup)
