import os
import interfaces/iSequence
import storage/containers/positionData
import storage/slidingDeque
import algorithm
import hts
import recordFilter
import jsonCollector
import consoleWriter
import chromosomeInjector

var bam: Bam
open(bam, paramStr(1), index=true)

var fai: Fai
if not open(fai, paramStr(2)):
  quit("Could not open fasta file.")

for chromosome in targets(bam.hdr):
  let name = chromosome.name

  var records = newRecordFilter(bam, name)
  var reference = fai.getISequence(name)
  var handler = newJsonCollector(newChromosomeInjector(newStreamWriter(), name))
  var storage = newSlidingDeque(200, handler.getICollector)

  pileup(records, reference, storage)
