## The module implements a storage mechanism for the ':pileup' data.  It uses a
## queue to organize 'PositionData' slots and collects data in a scanning
## left-to-right fasion. The storage can be flushed gradually.  This means that
## collected data can be submitted for further processing as soon as it finds
## itself behind the 'scan line', there is no need to wait for the whole pileup
## to finish. The module also provides a '%' (toJson) procedure which converts
## a 'PositionData' object into the standard LoFreq3 JSON format.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

import deques
import containers/positionData
import math
import ../pipetools
import ../../region


## Defines a type of the expected submit procedure. It should consume a
## 'PositionData' object without returning a result (since the storage should
## not be responsible for processing any data after it is collected.
type DataToVoid* =  proc(data: PositionData): void

## Defines a type which can be wrapped into a submit procedure type
## 'DataToVoid' with a discard operation.
type DataToType*[T] =  proc(data: PositionData): T


# NOTE: The queue does not really need to be double ended, but the 'Queue'
# module is deprecated.
type SlidingDeque* = ref object
  ## Defines a 'SlidingDeque' type and its relevant fields.
  deq: Deque[PositionData]
  submit: DataToVoid
  initialSize: int # estimated maximum size of the double ended queue
  beginning: int
  chromosome: string
  # Having the chromosome as a a field on the storage object is certainly less
  # than ideal. I will probably change this to be injected later.
  region: Region# FIXME this is a stupid hack to avoid submission of positions outside of region


const DEFAULT_INITIAL_SIZE = 200# FIXME autoset from readlength?


proc posWithinRegion(pos: PositionData, reg: Region): bool =
  if pos.refIndex <= reg.s or pos.refIndex > reg.e:
    return false
  else:
    return true


proc newSlidingDeque*(chromosome: string, region: Region, submit: DataToVoid,
                      initialSize: int = DEFAULT_INITIAL_SIZE
                     ): SlidingDeque {.inline.} =
  ## Constructs a new 'SlidingDeque' object. 
  ## The paramater 'submit' is a procedure expected to perform all furhter 
  ## processing. This procedure must be a consumer (not returning anything) of
  ## type 'DataToVoid'. If the user passes a mapping procedure of type 
  ## 'DataToType', it matches against the second constructor which performs the
  ## required wrapping.
  ## There is an optional initial size argument for the queue for optimization
  ## purposes.
  let adjustedSize = nextPowerOfTwo(initialSize)
  SlidingDeque(
    deq: initDeque[PositionData](adjustedSize),
    submit: submit,
    initialSize: adjustedSize,
    beginning: 0,
    chromosome: chromosome,
    region: region
  )


proc newSlidingDeque*(chromosome: string, submit: DataToType,
                      initialSize: int = DEFAULT_INITIAL_SIZE
                     ): SlidingDeque {.inline.} =
  ## Constructs a new 'SlidingDeque' object. 
  ## The paramater 'submit' is a procedure expected to perform all furhter 
  ## processing.
  ## This constructor accepts a procedure mapping a 'PositionData' object to
  ## any desired type and wraps it into a 'DataToVoid' consumer. There is an
  ## optional initial size argument for the queue for optimization purposes.
  newSlidingDeque(initialSize, chromosome, submit.done())


proc submitDeq(self: SlidingDeque,
               deq: var Deque[PositionData]): void {.inline.} =
  # FIXME: implement asyncronous procedure (probably outside of this module)
  # Submits the current deque for further processing
  for pd in deq:
    if posWithinRegion(pd, self.region):
      self.submit(pd)


proc resetDeq(self: SlidingDeque, beginning: int): void {.inline.} =
  # Submits all elements from the current deque for furhter processing.
  self.submitDeq(self.deq)
  self.deq = initDeque[PositionData](self.initialSize)
  self.beginning = beginning


proc `[]`(self: SlidingDeque, position:int): PositionData {.inline.} =
  ## Access a position in the deque, for testing purposes.
  if position < self.beginning or position >= self.beginning + self.deq.len:
    raise newException(ValueError, "Illegal position")
  return self.deq[position - self.beginning]


proc sanityCheck(beginning, length, position: int) : void {.inline.} =
  ## Checks whether the given position is valid based on the
  ## current deque beginning and length. Assumes that the deque is allowed
  ## to be extended.
  assert position >= beginning, "The file is not sorted: " &
    $position & ' ' & $beginning
  assert position <= (beginning + length), "Invalid position" &
    $position & $(beginning + length)


proc sanityCheckNoExtend(beginning, length, position: int): void {.inline.} =
  ## Checks whether the given position is valid based on the
  ## current deque beginning and length. Assumes that the deque is not allowed
  ## to be extended.
  assert position >= beginning, "The file is not sorted: " &
    $position & ' ' & $beginning
  assert position < beginning + length


proc ensureStorage(self: SlidingDeque, position:int,
                   refBase: char): void {.inline.} =
  ## Performs sanity checks before and, if needed, extends the storage.
  let length = self.deq.len
  sanityCheck(self.beginning, length, position)

  if position == (self.beginning + length):
    self.deq.addLast(newPositionData(position+1, refBase,
                        self.chromosome))


proc recordMatch*(self: SlidingDeque, position: int,
                  base: char, quality: int, reversed: bool,
                  refBase: char): void {.inline.} =
  ## Records match event information on for a given position.
  self.ensureStorage(position, refBase)
  self.deq[position - self.beginning].addMatch(base, quality, reversed)


proc recordDeletion*(self: SlidingDeque, position: int, bases: string,
                     quality: int, reversed: bool): void {.inline.} =
  ## Records deletion event information for a given position. If using this
  ## storage, all deletions should be reported on the base to their left. Thus,
  ## this procedure does not allow the deque to be extended and assumes the
  ## needed slot is already available.
  sanityCheckNoExtend(self.beginning, self.deq.len, position)
  self.deq[position - self.beginning].addDeletion(bases, quality, reversed)


proc recordInsertion*(self: SlidingDeque, position: int, bases: string,
                      quality: int, reversed: bool): void {.inline.} =
  ## Records insertion event infromation for a given position.
  sanityCheckNoExtend(self.beginning, self.deq.len, position)
  self.deq[position - self.beginning].addInsertion(bases, quality, reversed)


proc flushAll*(self: SlidingDeque): int {.inline.} =
  ## Submits all elements currently contained in the queue
  ## for further processing. The method returns the number of
  ## submitted elements.
  result = self.deq.len
  self.resetDeq(0)


proc flushUpTo*(self: SlidingDeque, position: int): int {.inline.} =
  ## Submits all elements with positions on the reference smaller than the 
  ## given argument for further processing. Enables the queue to slide.
  ## The method returns the number of submitted elements.
  if position < self.beginning - 1:
    raise newException(ValueError, "Flush index lower than beginning.")
  
  # if a new start position is larger than all positions contained in
  # the deque, instead of emptying it manually, we can submit it and
  # make a new one  
  if position >= self.beginning + self.deq.len:
    # FIXME: simplify once the
    # tests are stable, this is an overkill
    result = self.deq.len
    self.resetDeq(position)
    return result

  while self.beginning < position:
    let pd = self.deq.popFirst()
    if posWithinRegion(pd, self.region):
      self.submit(pd)
    self.beginning.inc
    result.inc


# when isMainModule:
  # block: # test constructor
  #   var pairs = [
  #     (given: 10, adjusted: 16),
  #     (given: 789, adjusted: 1024)
  #   ]
  #   for pair in pairs: 
  #     proc f(d: PositionData): void =
  #       discard
  #     var storage = newSlidingDeque(pair.given, f )
  #     doassert storage.submit == f
  #     doAssert storage.initialSize == pair.adjusted
  #     doAssert storage.deq.len == 0
  #     doAssert storage.beginning == 0

  #   block: # test record and get
  #     var acc: seq[PositionData] = @[]
  #     var storage = newSlidingDeque(10, proc (d: PositionData) = acc.add(d))
  #     storage.record(0, "A", 'A')
  #     storage.record(0, "C", 'A')
  #     storage.record(0, "C", 'A')
  #     storage.record(1, "G", 'G')
  #     storage.record(1, "-AA", 'T')
      
  #     doAssert storage.deq.len == 2
  #     doAssert storage.beginning == 0

  #     # legal gets
  #     let entry0 = storage[0]
  #     doAssert entry0.events["A"] == 1
  #     doAssert entry0.events["C"] == 2
  #     doAssert entry0.refBase == 'A'
  #     doAssert entry0.refIndex == 0

  #     let entry1 = storage[1]
  #     doAssert entry1.events["G"] == 1
  #     doAssert entry1.events["-AA"] == 1
  #     doAssert entry1.refBase == 'G'
  #     doAssert entry1.refIndex == 1

  #     # illegal gets
  #     try: 
  #       discard storage[2]
  #       doAssert false
  #     except ValueError:
  #       doAssert true
  #     except:
  #       doAssert false

  #     try: 
  #       discard storage[-1]
  #       doAssert false
  #     except ValueError:
  #       doAssert true
  #     except:
  #       doAssert false

  #   block: # test flush with overlap
  #     var acc: seq[PositionData] = @[]
  #     var storage = newSlidingDeque(20, proc (d: PositionData) = acc.add(d))
      
  #     # test legal record
  #     storage.recordMatch(0, 'A', 'A')
  #     storage.recordMatch(1, 'G', 'G')
  #     storage.recordMatch(2, 'A', 'T')
  #     doAssert storage.deq.len == 3
      
  #     doAssert storage.flushUpTo(2) == 2
      
  #     doAssert storage.beginning == 2
  #     doAssert storage.deq.len == 1
  #     doAssert acc.len == 2

  #     doAssert acc[0].refIndex == 0
  #     doAssert acc[0].refBase == 'A'

  #     doAssert acc[1].refBase == 'G'
  #     doAssert acc[1].refIndex == 1

  #   block: # test flush without overlap
  #     var acc: seq[PositionData] = @[]
  #     var storage = newSlidingDeque(20, proc (d: PositionData) = acc.add(d))
  #     storage.record(0, "A", 'A')
  #     storage.record(1, "G", 'G')
  #     storage.record(2, "A", 'T')

  #     doAssert storage.flushUpTo(4) == 3

  #     doAssert storage.beginning == 4
  #     doAssert storage.deq.len == 0
  #     doAssert acc.len == 3

  #     doAssert acc[0].refIndex == 0
  #     doAssert acc[0].refBase == 'A'

  #     doAssert acc[1].refBase == 'G'
  #     doAssert acc[1].refIndex == 1

  #     doAssert acc[2].refBase == 'T'
  #     doAssert acc[2].refIndex == 2

  #   block:
  #     var actual : seq[PositionData] = @[]
  #     var storage = newSlidingDeque(20, proc (d: PositionData): void = actual.add(d))
      
  #     assert storage.flushUpTo(0) == 0
  #     storage.record(0,"A", 'A')
  #     storage.record(1,"A", 'A')
  #     storage.record(2, "C", 'A')
  #     storage.record(3,"A", 'A')
  #     doAssert storage.deq.len == 4
  #     doAssert storage.beginning == 0

  #     assert storage.flushUpTo(0) == 0
  #     storage.record(0,"A", 'A')
  #     storage.record(1,"T", 'A')
  #     storage.record(2, "G", 'A')
  #     storage.record(3,"A", 'A')
  #     doAssert storage.deq.len == 4
  #     doAssert storage.beginning == 0

  #     assert storage.flushUpTo(0) == 0
  #     storage.record(0,"A", 'A')
  #     storage.record(1,"T", 'A')
  #     storage.record(2, "-AC", 'A')
  #     storage.record(3,"G", 'G')
  #     storage.record(4,"C", 'G')
  #     doAssert storage.deq.len == 5, $storage.deq.len
  #     doAssert storage.beginning == 0

  #     assert storage.flushUpTo(10) == 5
  #     storage.record(10, "A", 'A')
  #     doAssert storage.deq.len == 1
  #     doAssert storage.beginning == 10

  #     var expected = @[
  #       {"A": 3}.toTable, 
  #       {"A": 1, "T": 2}.toTable, 
  #       {"C": 1, "G": 1, "-AC": 1}.toTable,
  #       {"A": 2, "G": 1}.toTable,
  #       {"C": 1}.toTable
  #     ]

  #     for idx, pair in zip(actual, expected):
  #       echo pair[0][]
  #       echo pair[1]
  #       doAssert pair[0].events == pair[1], $pair[0].events & " " & $pair[1]

