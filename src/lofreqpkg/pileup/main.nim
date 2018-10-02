import os
import interfaces/iSequence
import storage/slidingDeque
import processor
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
  
  var injectChromosome = getJsonPropertyInjector("chromosome", name)
  var handler = toJson
    .thenDo(injectChromosome)
    .then(print)

  var storage = newSlidingDeque(200, handler)
  var processor = newProcessor(storage)
  pileup(records, reference, processor)
