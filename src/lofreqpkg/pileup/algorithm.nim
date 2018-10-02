## The module implements a pileup algorithm performed over a single chromosome and
## its reads.
import hts

proc allowed(operation: CigarOp): bool {.inline.} =
  ## Determines whether a cigar operation is allowed or not.
  ## Returns true if allowed, false otherwise.
  case operation
    of hard_clip, ref_skip, pad: false
    else: true


proc processEvent[TSequence, TProcessor](event: CigarElement, processor: var TProcessor,
                  read: Record, reference: TSequence,
                  readOffset: int, refOffset: int): (int, int) {.inline.} =
  ## Processes one event (cigar element) on the read and returns the updated offset 
  ## The event is described by the first argument, 'event'. The parameter 'processor'
  ## provides a reference to the processor used to record the events. 'read' and 'reference' are the
  ## read and the reference sequences. 'readOffset' and 'refOffset' mark the current position 
  ## on the read and the reference respectively. After processing the event, the function
  ## returns new offsets as a tuple.
  let operation = event.op
  
  if not operation.allowed:
    raise newException(ValueError, "Invalid operation: " & $operation)

  if operation == soft_clip:
    # a soft clip is not reported but advaces the position on the read
    return (readOffset + event.len, refOffset)
  
  let consumes = event.consumes()
  
  if consumes.query and consumes.reference:
    # mutual, report all matches
    processor.reportMatches(readOffset, refOffset, event.len,
                  read, reference)
    return (readOffset + event.len, refOffset + event.len)

  if consumes.reference:
    # reference only, report deletion
    processor.reportDeletion(readOffset, refOffset, event.len,
                   read, reference)
    return (readOffset, refOffset + event.len)

  if consumes.query:
    # read only, report insertion
    processor.reportInsertion(readOffset, refOffset, event.len,
                    read, reference)
    return (readOffset + event.len, refOffset)
  
  raise newException(ValueError, "Operation does not consume anything.")

proc valid(cigar: Cigar): bool {.inline.} =
  ## Tests whether a CIGAR is valid.
  ## Returns true if the CIGAR is valid, false otherwise.
  ## todo check the rest of the cigar to prevent later exceptions
  case cigar[0].op
    of CigarOp.insert, CigarOp.deletion: false
    else: true 

    
proc pileup*[TSequence, TRecordIterable, TProcessor](reads: TRecordIterable,
                                                 reference: TSequence,
                                                 processor: var TProcessor): void {.inline.} =
  ## Performs a pileup over all reads provided by the read iterable parameter and reports
  ## the results to the given processor object.
  ## The parameter `reads` is an iterable (has an `items` method) yielding all
  ## reads which should be piled up. It must yield objects of type 'Bam.Record'.
  ## The parameter `reference` is an interface for accessing the reference sequence.
  ## The parameter `processor` is an implementation of a processor object used in
  ## the pileup.
  for read in reads:
      let cigar = read.cigar

      if not cigar.valid:
        # Skipping all invalid reads (e.g. beginning with deletions/insertions)
        writeLine(stderr, "WARNING: Skipping read with invalid CIGAR: " & $read)
        continue

      var 
        readOffset = 0
        refOffset = read.start

      # since the file is sorted and a read CANNOT begin with an insertion or deletion,
      # we can safley flush any information related to
      # indices smaller than the current start of the read
      processor.beginRead(read.start)
      for event in cigar:
        (readOffset, refOffset) =
          processEvent(event, processor, read, reference,
                       readOffset, refOffset)

  # the pileup is done and all positions can be flushed      
  processor.done()


