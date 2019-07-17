## This module provides an implementation of the a processor used in the pileup
## algorithm (the module 'algorithm'). It gives the user an option to define
## provide functions for calculating base qualities of matches, insertions and
## deletions respecively, depending on specific needs.  The main idea behind
## separating the processing from the algorithm (iterating) is extensibility.
## This way, the 'Processor' can be configured to calculate qualities of all
## operation in arbitrary ways. Additionaly, if a user wants to substantially
## change the logic by which the data is processed (e.g. for a simpler form of
## pileup), they can implement a new processor and pass that to the algorithm
## instead without any need to change the algorithm itself.
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


## The expected type of procedures for calculating the qualities. All relevant
## information should be obtainable through the record and the index.
## NOTE: When dealing with a deletion or an insertion, the index is actually
## the index of the base TO THE LEFT (the reference index of the last mutual
## base.  Hovewer, this should not really make a difference since, in most
## cases, all deletion and insertion quality data is found in the record's
## custom fields.
type TQualityProc = proc (r: Record, i: int): int

type Processor[TStorage] = ref object
  ## The 'Processor' type. Its fields are configuration options.
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
  ## provided functions. The processor stores and updates the data in the
  ## provided storage object. All further handling is done by the storage
  ## object.
  Processor[TStorage](storage: storage,
                      matchQualityAt: matchQuality,
                      insertionQualityAt: insertionQuality,
                      deletionQualityAt: deletionQuality,
                     )


proc newProcessor*[TStorage](storage: TStorage):
                              Processor[TStorage] {.inline.} =
  ## The default constructor for the 'Processor' type.
  #FIXME use varargs to pass parameters down?
  Processor[TStorage](storage: storage,
                      matchQualityAt: proc(r: Record, i: int): int =
                        int(r.baseQualityAt(i)),
                      insertionQualityAt: proc(r: Record, i: int): int =
                        DEFAULT_INSERTION_QUALITY,
                      deletionQualityAt: proc(r: Record, i: int): int =
                        DEFAULT_DELETION_QUALITY)


proc processMatches*[TSequence](self: Processor,
                   readStart: int, refStart: int, length: int,
                   read: Record, reference: TSequence) : void {.inline.} =
  ## Processes a matching substring between the read and the reference. All
  ## necessary information is available through the arguments. A matching
  ## substring consists of multiple contiguous matching bases.
  for offset in countUp(0, length - 1):
    let refOff = refStart + offset
    let readOff = readStart + offset
    # q == -1 is filtering signal
    let bq = self.matchQualityAt(read, readOff)
    if bq != -1:
      self.storage.recordMatch(refOff, read.baseAt(readOff),
                               bq,
                               read.flag.reverse,
                               reference.baseAt(refOff))


proc processInsertion*[TSequence](self: Processor,
                     readStart: int, refIndex: int, length: int,
                     read: Record, reference: TSequence): void {.inline.} =
  ## Processes an insertion on the read (wrt. the reference). All necessary
  ## information is available through the arguments. An insertion consists of
  ## one or more bases found on the read, but not on the reference.
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
  ## Processes an deletion on the read (wrt. the reference). All necessary
  ## information is available through the arguments. A deletion consists of one
  ## or more bases found on the read, but not on the reference.
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
