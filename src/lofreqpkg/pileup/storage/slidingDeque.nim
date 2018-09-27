## The module implements a storage mechanism for the :pileup data.
## It uses a queue to organize 'PositionData' slots and it collects data
## in a scanning left-to-right fasion. The storage can be flushed gradually. 
## This means that data can be submitted for further processing as soon as it 
## finds itself behind the 'scan line', no need to wait for the whole pileup 
## to finish.
import deques
import containers/eventData
import containers/positionData
import math
import tables
import ../interfaces/pSubmit

# NOTE: The quue does not really need to be double ended, but the 'Queue' module
# is deprecated.
type SlidingDeque* = ref object
  deq: Deque[PositionData]
  submit: PSubmit[PositionData]
  initialSize: int # estimated maximum size of the double ended queue
  beginning: int


proc newSlidingDeque*(initialSize: int, submit: PSubmit[PositionData]): SlidingDeque =
  ## Constructs a new SlidingDeque object. 
  ## The paramater 'submit' is a function expected to perform all furhter processing.
  ## There is an optional initial size argument for the queue.
  let adjustedSize = nextPowerOfTwo(initialSize)
  SlidingDeque(
    deq: initDeque[PositionData](adjustedSize),
    submit: submit,
    initialSize: adjustedSize,
    beginning: 0
  )


proc submitDeq(self: SlidingDeque, deq: var Deque[PositionData]): void =
  # todo implement
  # asyncronous function (probably outside of this module)
  # Submits the current deque to another thread for furhter processing
  for element in deq:
    self.submit(element)


proc resetDeq(self: SlidingDeque, beginning: int) =
  # Submits all elements from the current deque for furhter processing.
  self.submitDeq(self.deq)
  self.deq = initDeque[PositionData](self.initialSize)
  self.beginning = beginning


proc `[]`(self: SlidingDeque, position:int): PositionData =
  ## Access a position in the deque, for testing purposes.
  if position < self.beginning or position >= self.beginning + self.deq.len:
    raise newException(ValueError, "Illegal position")
  return self.deq[position - self.beginning]


proc sanityCheck(beginning, length, position: int) : void =
  ## Checks whether the given position is valid based on the
  ## current deque beginning and length. Allows the deque to be extended. 
  assert position >= beginning, "The file is not sorted: " & $position & ' ' & $beginning
  assert position <= (beginning + length), "Invalid position" & $position & $(beginning + length)
  

proc sanityCheckNoExtend(beginning, length, position: int): void =
  ## Checks whether the given position is valid based on the
  ## current deque beginning and length. Does not allow the deque to be extended.
  assert position >= beginning, "The file is not sorted: " & $position & ' ' & $beginning
  assert position < beginning + length


proc extendStorage(self: SlidingDeque, position:int, refBase: char): void =
  let length = self.deq.len
  sanityCheck(self.beginning, length, position)

  if position == (self.beginning + length):
    self.deq.addLast(newPositionData(length + self.beginning, refBase))
  

proc recordMatch*(self: SlidingDeque, position: int,
                  base: char, quality: int, refBase: char): void =
  ## Records match event information on for a given position
  self.extendStorage(position, refBase)
  self.deq[position - self.beginning].addMatch(base, quality)
    

proc recordDeletion*(self: SlidingDeque, position: int,
                     bases: string, quality: int): void =
  ## Records deletion event information for a given position.
  sanityCheckNoExtend(self.beginning, self.deq.len, position)
  self.deq[position - self.beginning].addDeletion(bases, quality)


proc recordInsertion*(self: SlidingDeque, position: int,
                      bases: string, quality: int): void =
  ## Records insertion event infromation for a given position.
  sanityCheckNoExtend(self.beginning, self.deq.len, position)
  self.deq[position - self.beginning].addInsertion(bases, quality)


proc flushAll*(self: SlidingDeque): int =
  ## Submits all elements currently contained in the queue 
  ## for further processing. The method returns the number of
  ## submitted elements
  result = self.deq.len
  self.resetDeq(0)


proc flushUpTo*(self: SlidingDeque, position: int): int =
  ## Submits all elements with positions on the reference smaller than the given
  ## argument for further processing. Enables the queue to slide. The method returns number
  ## of submitted elements.
  if position < self.beginning - 1:
    raise newException(ValueError, "Flush index lower than beginning.")
  
  # if a new start position is larger than all positions contained in
  # the deque, instead of emptying it manually, we can submit it and
  # make a new one  
  if position >= self.beginning + self.deq.len:
    result = self.deq.len
    self.resetDeq(position)
    self.beginning = position
    return result

  while self.beginning < position:
    self.submit(self.deq.popFirst())
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
  #     doAssert entry0.referenceBase == 'A'
  #     doAssert entry0.referenceIndex == 0

  #     let entry1 = storage[1]
  #     doAssert entry1.events["G"] == 1
  #     doAssert entry1.events["-AA"] == 1
  #     doAssert entry1.referenceBase == 'G'
  #     doAssert entry1.referenceIndex == 1

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

  #     doAssert acc[0].referenceIndex == 0
  #     doAssert acc[0].referenceBase == 'A'

  #     doAssert acc[1].referenceBase == 'G'
  #     doAssert acc[1].referenceIndex == 1

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

  #     doAssert acc[0].referenceIndex == 0
  #     doAssert acc[0].referenceBase == 'A'

  #     doAssert acc[1].referenceBase == 'G'
  #     doAssert acc[1].referenceIndex == 1

  #     doAssert acc[2].referenceBase == 'T'
  #     doAssert acc[2].referenceIndex == 2

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

