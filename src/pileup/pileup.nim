## Implementation of an alternative procedure to samtools' pileup that does not require 
## keeping pointers to all the reads in memory and provides support for streaming results
## to microservices as soon as possible. This was inspired by https://brentp.github.io/post/no-pile/.
import hts 
import math

# I should perhaps remove this as it is just a special case for 
# reportMatches and one additional function call per every match may
# drastically reduce performance. However, I left it for semantic reasons and encapsulation
# (e.g. beginning of the read). 
# todo refactor after benchmarking
proc reportMatch[TSequence, TStorage](storage: var TStorage,
                 readIndex: int, refIndex: int,
                 read: Record, reference: TSequence): void =
  ## Reports one matching base between the read and the reference to the given
  ## storage object.
  storage.record(refIndex, $read.baseAt(readIndex), reference.baseAt(refIndex))



proc reportMatches[TSequence, TStorage](storage: var TStorage, 
                   readStart: int, refStart: int, length: int, 
                   read: Record, reference: TSequence) : void =
  ## Reports a matching substring between the read and the reference to
  ## the given storage object.
  ## A matching substring consists of multiple continuos matching bases.
  if length == 1:
    reportMatch(storage, readStart, refStart, read, reference)
    return

  for offset in countUp(0, length - 1):
    reportMatch(storage, readStart + offset, refStart + offset, read, reference)
  

proc reportInsertion[TSequence, TStorage](storage: var TStorage,
                     readStart: int, refIndex: int, length: int,
                     read: Record, reference: TSequence): void =
  ## Reports an insertion on the read (wrt. the reference) to the provided storage.
  ## An insertion consists of one or more bases found on the read,
  ## but not on the reference.
  var value = "+"

  for offset in countUp(readStart, readStart + length - 1):
    value &= read.baseAt(offset)

  # insertion is reported on the base that preceeds it
  storage.record(refIndex - 1, value , '/')



proc reportDeletion[TSequence, TStorage](storage: var TStorage,
                    readStart: int, refStart: int, length: int,
                    reference: TSequence): void =
  ## Reports a deletion on the read (wrt. the reference) to the provided storage.\
  ## A deletion consists of one or more bases found on the reference,
  ## but not on the reads.
  var value = "-"

  for offset in countUp(refStart, refStart + length - 1):
    value &= reference.baseAt(offset) 
    storage.record(offset, "*", reference.baseAt(offset))

  # deletion is reported on the base that preceeds it
  storage.record(refStart - 1, value, '/')



proc processEvent[TSequence, TStorage](event: CigarElement, storage: var TStorage, 
                  read: Record, reference: TSequence,
                  readOffset: var int, refOffset: var int): void =
  ## Processes one event (cigar element) on the read

  let operation = event.op
  if operation == soft_clip:
    readOffset += event.len
    return
  if operation == hard_clip: raise newException(ValueError, "hard clip")
  assert operation != ref_skip and operation != pad, "Illegal operation"


  let consumes = event.consumes()
  
  if consumes.query and consumes.reference: 
    # mutual, report all matches
    reportMatches(storage, 
                  readOffset, refOffset, event.len, 
                  read, reference)
    readOffset += event.len
    refOffset += event.len

  elif consumes.reference:
    # reference only, report deletion
    reportDeletion(storage, 
                   readOffset, refOffset, event.len, 
                   reference)
    refOffset += event.len

  elif consumes.query:
    # read only, report insertion
    reportInsertion(storage, 
                    readOffset, refOffset, event.len,
                    read, reference)
    readOffset += event.len

proc pileup*[TSequence, TReadIterable, TStorage](reads: TReadIterable, 
                                                 reference: TSequence,
                                                 storage: var TStorage): void =
  ## Performs a pileup over all reads provided by the read iterable parameter and reports
  ## the results to the given storage object.
  ## The parameter `reads` is an iterable (has an `items` method) providing all
  ## reads which should be piled up.
  ## The parameter `reference` is an interface for accessing the reference sequence.
  ## The parameter `storage` is an implementation of a storage object used in
  ## the pileup.
  for read in reads:
      var 
        readOffset = 0
        refOffset = read.start

      # since the file is sorted and a read CANNOT begin with an insertion or deletion,
      # we can safley flush any information related to
      # indices smaller than the current start of the read
      var a = storage.flushUpTo(read.start- 1)
      for event in read.cigar:
        processEvent(event, storage, read, reference, 
                     readOffset, refOffset)
  discard storage.flushAll()
