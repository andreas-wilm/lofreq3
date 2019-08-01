## The module implements a pileup algorithm performed over a single chromosome
## and its reads. It serves as a template and defines only the basic procedure
## for iterating over the CIGAR strings of BAM records. The processing itself
## is defined by an outside object.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License


# standard
import os
import logging
# third party
import hts
# project specific
import ../region
import recordFilter
import interfaces/iSequence
import storage/slidingDeque
import processor
import storage/slidingDeque


var logger = newConsoleLogger(fmtStr = verboseFmtStr, useStderr = true)

        
proc allowed(operation: CigarOp): bool {.inline.} =
  ## Determines whether a cigar operation is allowed or not.
  ## Returns true if allowed, false otherwise.
  case operation
    of hard_clip, ref_skip, pad: false
    else: true


proc processEvent[TSequence, TProcessor](event: CigarElement,
                  processor: var TProcessor,
                  read: Record, reference: TSequence,
                  readOffset: int, refOffset: int): (int, int) {.inline.} =
  ## Processes one event (cigar element) on the read and returns the updated
  ## offset. The event is described by the first argument, 'event'.
  ## The parameter 'processor' provides a reference to the processor used to
  ## handle the events. 'read' and 'reference' are the read and the reference
  ## sequences. 'readOffset' and 'refOffset' mark the current position on the 
  ## read and the reference respectively. After processing the event, the 
  ## procedure returns the new offsets as a tuple.
  let operation = event.op

  if not operation.allowed:
    raise newException(ValueError, "Invalid operation: " & $operation)

  if operation == soft_clip:
    # a soft clip is not processed but advaces the position on the read
    return (readOffset + event.len, refOffset)

  let consumes = event.consumes()

  if consumes.query and consumes.reference:
    # mutual, process all matches
    processor.processMatches(readOffset, refOffset, event.len,
                    read, reference)
    return (readOffset + event.len, refOffset + event.len)

  if consumes.reference:
    # reference only, process deletion
    processor.processDeletion(readOffset, refOffset, event.len,
                     read, reference)
    return (readOffset, refOffset + event.len)

  if consumes.query:
    # read only, process insertion
    processor.processInsertion(readOffset, refOffset, event.len,
                      read, reference)
    return (readOffset + event.len, refOffset)

  raise newException(ValueError, "Operation does not consume anything.")


proc valid(cigar: Cigar): bool {.inline.} =
  ## Tests whether a CIGAR is valid.
  ## Returns true if the CIGAR is valid, false otherwise.
  ## todo check the rest of the cigar to prevent later exceptions
  result = true
  if cigar[0].op in [CigarOp.insert, CigarOp.deletion]:
    result = false
  elif cigar[0].op == CigarOp.soft_clip:
    if len(cigar)>1 and cigar[1].op in [CigarOp.insert, CigarOp.deletion]:
        result = false


proc pileup*(fai: Fai, records: RecordFilter, region: Region, handler: DataToVoid): void {.inline.} =
  ## Performs a pileup over all reads provided by records

  var reference: ISequence# our own type, hence using loadSequence below
  var storage = newSlidingDeque(records.chromosomeName, region, handler)
  var processor = newProcessor(storage)

  for read in records:
      # all records come from the same chromosome as guarantted by RecordFilter
      # load reference only after we're sure there's data to process
      if reference.len == 0:
        reference = fai.loadSequence(records.chromosomeName)
       
      let cigar = read.cigar

      if not cigar.valid:
        # Skipping all invalid reads
        logger.log(lvlWarn, "Skipping read with invalid CIGAR: " & $read)
        continue

      var
        readOffset = 0
        refOffset = read.start

      # tell the processor that the new read is about to start
      processor.beginRead(read.start)

      # process all events on the read
      for event in cigar:
        (readOffset, refOffset) =
          processEvent(event, processor, read, reference,
                       readOffset, refOffset)

  # inform the processor that the pileup is done
  processor.done()

