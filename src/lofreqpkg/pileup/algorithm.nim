## The module implements a pileup algorithm performed over a single chromosome and
## its reads.
import hts

# placeholders until we find a way to record true qualities
const DEFAULT_DELETION_QUALITY = 40
const DEFAULT_INSERTION_QUALITY = 40
const DEFAULT_BLANK_QUALITY = 40
const DEFAULT_BLANK_SYMBOL = '*'

proc allowed(operation: CigarOp): bool =
  ## Determines whether a cigar operation is allowed or not.
  ## Returns true if allowed, false otherwise.
  case operation
    of hard_clip, ref_skip, pad: false
    else: true


proc reportMatches[TSequence, TStorage](storage: var TStorage,
                   readStart: int, refStart: int, length: int,
                   read: Record, reference: TSequence) : void =
  ## Reports a matching substring between the read and the reference to
  ## the given storage object.
  ## A matching substring consists of multiple continuos matching bases.
  for offset in countUp(0, length - 1):
    let refOff = refStart + offset
    let readOff = readStart + offset
    storage.recordMatch(refOff, read.baseAt(readOff), int(read.baseQualityAt(readOff)),
                        reference.baseAt(refOff))

  

proc reportInsertion[TSequence, TStorage](storage: var TStorage,
                     readStart: int, refIndex: int, length: int,
                     read: Record, reference: TSequence): void =
  ## Reports an insertion on the read (wrt. the reference) to the provided storage.
  ## An insertion consists of one or more bases found on the read,
  ## but not on the reference.
  var value = ""
  for offset in countUp(readStart, readStart + length - 1):
    value &= read.baseAt(offset)

  # insertion is reported on the base that preceeds it
  storage.recordInsertion(refIndex - 1, value, DEFAULT_INSERTION_QUALITY)



proc reportDeletion[TSequence, TStorage](storage: var TStorage,
                    readStart: int, refStart: int, length: int,
                    read: Record, reference: TSequence): void =
  ## Reports a deletion on the read (wrt. the reference) to the provided storage.
  ## A deletion consists of one or more bases found on the reference,
  ## but not on the reads.
  var value = ""
  for offset in countUp(refStart, refStart + length - 1):
    value &= reference.baseAt(offset)
    storage.recordMatch(offset, DEFAULT_BLANK_SYMBOL ,DEFAULT_BLANK_QUALITY,
                        reference.baseAt(offset))

  # deletion is reported on the base that preceeds it
  storage.recordDeletion(refStart - 1, value, DEFAULT_DELETION_QUALITY)



proc processEvent[TSequence, TStorage](event: CigarElement, storage: var TStorage,
                  read: Record, reference: TSequence,
                  readOffset: int, refOffset: int): (int, int) =
  ## Processes one event (cigar element) on the read and returns the updated offset 
  ## The event is described by the first argument, 'event'. The parameter 'storage'
  ## provides a reference to the storage used to record the events. 'read' and 'reference' are the
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
    reportMatches(storage,
                  readOffset, refOffset, event.len,
                  read, reference)
    return (readOffset + event.len, refOffset + event.len)

  if consumes.reference:
    # reference only, report deletion
    reportDeletion(storage,
                   readOffset, refOffset, event.len,
                   read, reference)
    return (readOffset, refOffset + event.len)

  if consumes.query:
    # read only, report insertion
    reportInsertion(storage,
                    readOffset, refOffset, event.len,
                    read, reference)
    return (readOffset + event.len, refOffset)
  
  raise newException(ValueError, "Operation does not consume anything.")

proc valid(cigar: Cigar): bool =
  ## Tests whether a CIGAR is valid.
  ## Returns true if the CIGAR is valid, false otherwise.
  ## todo check the rest of the cigar to prevent later exceptions
  case cigar[0].op
    of CigarOp.insert, CigarOp.deletion: false
    else: true 

    
proc pileup*[TSequence, TRecordIterable, TStorage](reads: TRecordIterable,
                                                 reference: TSequence,
                                                 storage: var TStorage): void =
  ## Performs a pileup over all reads provided by the read iterable parameter and reports
  ## the results to the given storage object.
  ## The parameter `reads` is an iterable (has an `items` method) yielding all
  ## reads which should be piled up. It must yield objects of type 'Bam.Record'.
  ## The parameter `reference` is an interface for accessing the reference sequence.
  ## The parameter `storage` is an implementation of a storage object used in
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
      discard storage.flushUpTo(read.start)
      for event in cigar:
        (readOffset, refOffset) =
          processEvent(event, storage, read, reference,
                       readOffset, refOffset)

  # the pileup is done and all positions can be flushed      
  discard storage.flushAll()


