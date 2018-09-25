import os
import interfaces/iSequence
import storage/containers/positionData
import storage/slidingDeque
import messaging
import pileup
import hts

type RecordContainer = ref object
  bam: Bam
  chromosomeName: string

proc newRecordContainer(bam: Bam, chromosomeName: string): RecordContainer =
  return RecordContainer(bam: bam, chromosomeName: chromosomeName)

iterator items(self: RecordContainer) : Record = 
  for read in self.bam.querys(self.chromosomeName):
    yield read

var bam: Bam
open(bam, paramStr(1), index=true)

var fai: Fai
if not open(fai, paramStr(2)):
  quit("Could not open fasta file.")

for chromosome in targets(bam.hdr):
  var storage = newslidingDeque(
                  200,
                  proc(d: PositionData): void = writeLine(stdout, createJsonMessage(d)) 
                )
  var records = newRecordContainer(bam, chromosome.name)
  pileup(records, fai.getISequence(chromosome.name), storage)
