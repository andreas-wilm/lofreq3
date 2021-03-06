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

const
  DEFAULT_MIN_COV* = 1
  DEFAULT_MAX_COV* = high(int)
  DEFAULT_MIN_BQ* = 3
  DEFAULT_USE_MQ* = true

type PileupParams* = object
  minCov*: Natural
  maxCov*: Natural
  minBQ*: Natural
  useMQ*: bool
  # FIXME add regions to plpParams


var plpParams*: PileupParams
plpParams = PileupParams(minCov: DEFAULT_MIN_COV,
                         maxCov: DEFAULT_MAX_COV,
                         minBQ: DEFAULT_MIN_BQ,
                         useMQ: DEFAULT_USE_MQ)

var logger = newConsoleLogger(fmtStr = verboseFmtStr, useStderr = true)


proc allowed(operation: CigarOp): bool {.inline.} =
  ## Determines whether a cigar operation is allowed or not.
  ## Returns true if allowed, false otherwise.
  case operation
    of hard_clip, ref_skip, pad: false
    else: true


proc processEvent[TSequence, TProcessor](event: CigarElement,
                  nextevent: CigarElement,
                  processor: var TProcessor,
                  read: Record, reference: TSequence,
                  readOffset: int, refOffset: int64): (int, int64) {.inline.} =
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
    # a soft clip is not processed but advances the position on the read
    return (readOffset + event.len, refOffset)

  let consumes = event.consumes()

  if consumes.query and consumes.reference:
    # mutual, process all matches
    processor.processMatches(readOffset, refOffset, event.len,
      read, reference, nextevent)
    return (readOffset + event.len, refOffset + event.len)

  elif consumes.reference:# also true for 'N'
    # reference only, process deletion
    if event.op == CigarOp.deletion:
      processor.processDeletion(readOffset, refOffset, event.len,
        read, reference)
    return (readOffset, refOffset + event.len)

  elif consumes.query:# also true for 'S'
    # read only, process insertion
    if event.op == CigarOp.insert:
      processor.processInsertion(readOffset, refOffset, event.len,
        read, reference)
    return (readOffset + event.len, refOffset)

  else:
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


proc pileup*(fai: Fai, records: RecordFilter, region: Region,
             handler: DataToVoid): void {.inline.} =
  ## Performs a pileup over all reads provided by records

  var reference: ISequence# our own type, hence using loadSequence below
  var storage = newSlidingDeque(records.chromosomeName, region, handler,
    plpParams.mincov, plpParams.maxcov)
  var processor = newProcessor(storage, plpParams.useMQ, plpParams.minBQ)
  
  for read in records:
      let cigar = read.cigar
      if not cigar.valid:
        # Skipping all invalid reads
        logger.log(lvlWarn, "Skipping read with invalid CIGAR: " & $read)
        continue
      
      # all records come from the same chromosome as guaranteed by RecordFilter
      # load reference only after we're sure there's data to process
      if reference.len == 0:
        reference = fai.loadSequence(records.chromosomeName)

      var
        readOffset = 0
        refOffset = int64(read.start)
      
      processor.beginRead(read)

      # process all events on the read. unfortunately we need to know
      # the next event to avoid storing indel quals twice, which makes
      # this a bit ugly
      for idx in 0..<len(cigar):
        let event = cigar[idx]
        var nextevent: CigarElement
        # uninitialized nextevent (= last element) translates to 0M
        if idx+1 < len(cigar):
          nextevent = cigar[idx+1]
        (readOffset, refOffset) = processEvent(event, nextevent,
          processor, read, reference, readOffset, refOffset)

  # inform the processor that the pileup is done
  processor.done()

