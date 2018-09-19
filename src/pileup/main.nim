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

func lenAsInt(target: Target): int =
  result = cast[int](target.length)
  # todo check integer limits

var bam: Bam
open(bam, paramStr(1), index=true)

var fai: Fai
if not open(fai, paramStr(2)):
  quit("Could not open fasta file.")

for chromosome in targets(bam.hdr):
  var storage = newslidingDeque(
                  20,
                  proc(d: PositionData): void = echo createJsonMessage(d)
                )
  var records = newRecordContainer(bam, chromosome.name)
  pileup(records, fai.getISequence(), storage)
