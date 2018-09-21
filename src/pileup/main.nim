import os
import interfaces/iSequence
import storage/containers/positionData
import storage/slidingDeque
import messaging
import pileup
import hts
import recordFilter

var bam: Bam
open(bam, paramStr(1), index=true)

var fai: Fai
if not open(fai, paramStr(2)):
  quit("Could not open fasta file.")

for chromosome in targets(bam.hdr):
  var records = newRecordFilter(bam, chromosome.name)
  var reference = fai.getISequence(chromosome.name)
  var storage = newslidingDeque(
                  200,
                  proc(d: PositionData): void = writeLine(stdout, createJsonMessage(d)) 
                )
  pileup(records, reference, storage)
