import os
import interfaces/iSequence
import storage/slidingDeque
import algorithm
import hts
import recordFilter
import util
import postprocessing

var bam: Bam
open(bam, paramStr(1), index=true)

var fai: Fai
if not open(fai, paramStr(2)):
  quit("Could not open fasta file.")

for chromosome in targets(bam.hdr):
  let name = chromosome.name

  var records = newRecordFilter(bam, name)
  var reference = fai.loadSequence(name)
  
  var injectChromosome = chromosomeInjector(name)
  var handler = toJson
    .then(injectChromosome)
    .thenDo(printOutput)

  var storage = newSlidingDeque(200, handler)
  
  pileup(records, reference, storage)
