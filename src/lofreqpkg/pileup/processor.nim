## This module provides an implementation of the a processor used in the
## pileup algorithm (the module 'algorithm'). It gives the user an option
## to define provide functions for calculating base qualities of matches,
## insertions and deletions respecively, depending on specific needs. 
##
## - Author: Filip Sodic <filip.sodic@gmail.com>
## - License: The MIT License
import hts
import strutils

# placeholders until we find a way to record true qualities
const DEFAULT_DELETION_QUALITY = 40
const DEFAULT_INSERTION_QUALITY = 40
const DEFAULT_BLANK_QUALITY = -1 
const DEFAULT_BLANK_SYMBOL = '*' # missing position symbol

## The expected type of procedures for calculating the qualities. All
## relevant information should be obtainable through the record and the index.
## NOTE: When dealing with a deletion or an insertion, the index is actually 
## the index of the base TO THE LEFT (the reference index of the last mutual
## base. Hovewer, this should not really make a difference since, in most
## cases, all deletion and insertion quality data is found in the record's
## custom fields.
type TQualityProc = proc (r: Record, i: int): int 


## The 'Processor' type. Its fields are configuration options.
type Processor[TStorage] = ref object
  storage: TStorage
  matchQualityAt: TQualityProc
  deletionQualityAt: TQualityProc
  insertionQualityAt: TQualityProc


proc newProcessor*[TStorage](storage: TStorage,
                             matchQuality: TQualityProc,
                             insertionQuality: TQualityProc,
                             deletionQuality: TQualityProc
                            ): Processor[TStorage] {.inline.} =
  ## Creates a new 'Processor'. The operation quality is calculated using the
  ## provided functions. The the processor stores and updates the data to the 
  ## provided storage object. All further handling is done by the storage object.
  Processor[TStorage](storage: storage,
                      matchQualityAt: matchQuality,
                      insertionQualityAt: insertionQuality,
                      deletionQualityAt: deletionQuality
                     )

proc newProcessor*[TStorage](storage: TStorage):
                              Processor[TStorage] {.inline.} =
  ## The default constructor for the 'Processor' type.
  Processor[TStorage](storage: storage,
                      matchQualityAt: proc(r: Record, i: int): int =
                        int(r.baseQualityAt(i)),
                      insertionQualityAt: proc(r: Record, i: int): int =
                        DEFAULT_INSERTION_QUALITY,
                      deletionQualityAt: proc(r: Record, i: int): int =
                        DEFAULT_DELETION_QUALITY
                     )


proc processMatches*[TSequence](self: Processor,
                   readStart: int, refStart: int, length: int,
                   read: Record, reference: TSequence) : void {.inline.} =
  ## Reports a matching substring between the read and the reference to
  ## the given storage object.
  ## A matching substring consists of multiple continuos matching bases.
  for offset in countUp(0, length - 1):
    let refOff = refStart + offset
    let readOff = readStart + offset
    self.storage.recordMatch(refOff, read.baseAt(readOff),
                             self.matchQualityAt(read, readOff),
                             read.flag.reverse,
                             reference.baseAt(refOff))

  

proc processInsertion*[TSequence](self: Processor,
                     readStart: int, refIndex: int, length: int,
                     read: Record, reference: TSequence): void {.inline.} =
  ## Reports an insertion on the read (wrt. the reference) 
  ## to the provided storage.
  ## An insertion consists of one or more bases found on the read,
  ## but not on the reference.
  var value = ""
  for offset in countUp(readStart, readStart + length - 1):
    value &= read.baseAt(offset)

  # insertion is reported on the base that preceeds it
  self.storage.recordInsertion(refIndex - 1, value,
                               self.insertionQualityAt(read, readStart),
                               read.flag.reverse)



proc processDeletion*[TSequence](self: Processor,
                    readIndex: int, refStart: int, length: int,
                    read: Record, reference: TSequence): void {.inline.} =
  ## Reports a deletion on the read (wrt. the reference) to the 
  ## provided storage.
  ## A deletion consists of one or more bases found on the reference,
  ## but not on the reads.
  var value = ""
  for offset in countUp(refStart, refStart + length - 1):
    value &= reference.baseAt(offset)
    self.storage.recordMatch(offset, DEFAULT_BLANK_SYMBOL,
                             DEFAULT_BLANK_QUALITY,
                             read.flag.reverse,
                             reference.baseAt(offset))

  # deletion is reported on the base that preceeds it
  self.storage.recordDeletion(refStart - 1, value,
                              self.deletionQualityAt(read, readIndex),
                              read.flag.reverse)


proc beginRead*(self: Processor, start: int): void {.inline.} =
  ## Performs what is necessary before starting a new read. In this case,
  ## this means flushing the storage up to the starting position.
  discard self.storage.flushUpTo(start)

proc done*(self: Processor): void {.inline.} =
  ## Finishes the processing, flushes the entire storage.
  discard self.storage.flushAll()
