## Implementation of an alternative procedure to samtools' pileup that does not require 
## keeping pointers to all the reads in memory and provides support for streaming results
## to microservices as soon as possible. This was inspired by https://brentp.github.io/post/no-pile/.
import os
import hts
import interfaces/iSequence
import storage/slidingDeque
import recordFilter
import algorithm
import util
import postprocessing


proc pileup*(bamFname: string, faFname: string) =
  var bam: Bam
  var fai: Fai

  if not open(fai, faFname):
    quit("Could not open fasta file.")
  
  open(bam, bamFname, index=true)

  for chromosome in targets(bam.hdr):
    let name = chromosome.name

    var records = newRecordFilter(bam, name)
    var reference = fai.loadSequence(name)
    
    var injectChromosome = getJsonPropertyInjector("chromosome", name)
    var handler = toJson
      .thenDo(injectChromosome)
      .then(print)

    var storage = newSlidingDeque(200, handler)
    
    pileup(records, reference, storage)

when isMainModule:
  import cligen
  dispatch(pileup)
