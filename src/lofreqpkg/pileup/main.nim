import os
import hts
import interfaces/iSequence
import storage/slidingDeque
import processor
import recordFilter
import algorithm
import pipetools
import postprocessing
import ../utils.nim


proc sampleQualAt(r: Record, i: int): int =
  let be = qual2Prob(int(r.baseQualityAt(i)))
  let me = qual2Prob(int(r.mappingQuality()))

  let je = me + (1 - me) * be
  result = prob2Qual(je)
  

var bam: Bam
open(bam, paramStr(1), index=true)

var fai: Fai
if not open(fai, paramStr(2)):
  quit("Could not open fasta file.")

for chromosome in targets(bam.hdr):
  let name = chromosome.name

  var records = newRecordFilter(bam, name)
  var reference = fai.loadSequence(name)
  
  var handler = toJson
    .then(print)

  var storage = newSlidingDeque(name, handler)
  var processor = newProcessor(storage,
                               sampleQualAt,
                               proc(r:Record, i: int): int = 46,
                               proc(r:Record, i: int): int = 46)
  pileup(records, reference, processor)
