## The module implements a pileup algorithm performed over a single chromosome
## and its reads. It serves as a template and defines only the basic procedure
## for iterating over the CIGAR strings of BAM records. The processing itself
## is defined by an outside object. 
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License


import hts

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
  ## function returns the new offsets as a tuple.
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
  case cigar[0].op
    of CigarOp.insert, CigarOp.deletion: false
    else: true 

    
proc pileup*[TSequence, TRecordIterable, TProcessor](reads: TRecordIterable,
                                                     reference: TSequence,
                                                     processor: var TProcessor
                                                    ): void {.inline.} =
  ## Performs a pileup over all reads provided by the read iterable parameter 
  ## and the passes all relevant information to the provided 'processor' object
  ## which is in charge of processing the operations.
  ## The parameter `reads` is an iterable (has an `items` method) yielding all
  ## reads which to be piled up. It must yield objects of type 'Bam.Record'.
  ## The parameter `reference` is an interface for accessing the reference
  ## sequence.
  ## The parameter `processor` must provide an implementation of the event 
  ## processing methods. 
  for read in reads:
      let cigar = read.cigar

      if not cigar.valid:
        # Skipping all invalid reads (e.g. beginning with deletions/insertions)
        writeLine(stderr, "WARNING: Skipping read with invalid CIGAR: " & $read)
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

