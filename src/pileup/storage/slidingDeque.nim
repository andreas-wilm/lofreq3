import deques
import containers/eventData
import containers/positionData
import math
import sequtils
import ../interfaces/iStorage
import tables 

type SlidingDeque* = ref object
  deq: Deque[PositionData]
  submit: proc (data: PositionData): void
  initialSize: int # estimated maximum size of the double ended queue
  beginning: int

proc newSlidingDeque*(initialSize: int, 
                    submitProc: proc (data: PositionData): void): SlidingDeque =
  let adjustedSize = nextPowerOfTwo(initialSize)
  SlidingDeque(
    deq: initDeque[PositionData](adjustedSize),
    submit: submitProc, 
    initialSize: adjustedSize, 
    beginning: 0
  )

proc submitDeq(self: SlidingDeque, deq: var Deque[PositionData]): void =
  # todo implement
  # asyncronous function
  # Submits the current deque to another thread for handling
  for element in deq:
    self.submit(element)


proc resetDeq(self: SlidingDeque, beginning: int) =
  # Submits the current deque and replaces it with a new one
  #
  # PARAMTERS:
  # self - this sliding deque
  # beginning - the beginning position of the new deque
  self.submitDeq(self.deq)
  self.deq = initDeque[PositionData](self.initialSize)
  self.beginning = beginning

proc `[]`(self: SlidingDeque, position:int): PositionData =
  if position < self.beginning or position >= self.beginning + self.deq.len:
    raise newException(ValueError, "Illegal position")
  return self.deq[position - self.beginning]

proc record*(
    self: SlidingDeque, 
    position: int,
    value: string,
    refBase: char
  ): void =
  let length = self.deq.len

  assert position >= self.beginning, "The file is not sorted: " & $position & ' ' & $self.beginning
  assert position <= (self.beginning + length), "Invalid position" & $position & $(self.beginning + length)
  
  if position == (self.beginning + length):
    self.deq.addLast(newPositionData(length + self.beginning, refBase))
  
  self.deq[position - self.beginning].increment(value)

proc flushAll*(self: SlidingDeque): int =
  result = self.deq.len
  self.resetDeq(0)

proc flushUpTo*(self: SlidingDeque, position: int): int = 
  ## Flushes/submits all finished slots in the storage. This is meant to be 
  ## called when starting to process a new read
  ## 
  ## @return the number of flushed items (primarily for testing and debugging purposes,
  ## feel free to remove when in release)
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

  while self.beginning < position: # -1 because of the insertions on read starts
    self.submit(self.deq.popFirst())
    self.beginning.inc
    result.inc

proc recordStart(
    self: SlidingDeque, 
    position: int, 
    readValue: string,
    refBase: char
  ): void =
  ## Used to handle the first position in a new read.
  ## Submits all positions which all smaller and only then updates
  ## the storage. 
  ## 
  ## PARAMETERS:
  ## self - this sliding deque
  ## position - the starting index of the new read (wrt. the reference)
  ## readValue - the string found at the position on the read
  discard self.flushUpTo(position)
  self.record(position, readValue, refBase)

proc getIStorage*(self: SlidingDeque): IStorage=
  return (
        record: proc(position: int,value: string,refBase: char): void= 
          self.record(position, value, refBase),
        flushUpTo: proc(position: int): int = 
          self.flushUpTo(position),
        flushAll: proc(): int  = 
          self.flushAll()
       )

when isMainModule:
  block: # test constructor
    var pairs = [
      (given: 10, adjusted: 16),
      (given: 789, adjusted: 1024)
    ]
    for pair in pairs: 
      proc f(d: PositionData): void =
        discard
      var storage = newSlidingDeque(pair.given, f )
      doassert storage.submit == f
      doAssert storage.initialSize == pair.adjusted
      doAssert storage.deq.len == 0
      doAssert storage.beginning == 0

    block: # test record and get
      var acc: seq[PositionData] = @[]
      var storage = newSlidingDeque(10, proc (d: PositionData) = acc.add(d))
      storage.record(0, "A", 'A')
      storage.record(0, "C", 'A')
      storage.record(0, "C", 'A')
      storage.record(1, "G", 'G')
      storage.record(1, "-AA", 'T')
      
      doAssert storage.deq.len == 2
      doAssert storage.beginning == 0

      # legal gets
      let entry0 = storage[0]
      doAssert entry0.events["A"] == 1
      doAssert entry0.events["C"] == 2
      doAssert entry0.referenceBase == 'A'
      doAssert entry0.referenceIndex == 0

      let entry1 = storage[1]
      doAssert entry1.events["G"] == 1
      doAssert entry1.events["-AA"] == 1
      doAssert entry1.referenceBase == 'G'
      doAssert entry1.referenceIndex == 1

      # illegal gets
      try: 
        discard storage[2]
        doAssert false
      except ValueError:
        doAssert true
      except:
        doAssert false

      try: 
        discard storage[-1]
        doAssert false
      except ValueError:
        doAssert true
      except:
        doAssert false

    block: # test flush with overlap
      var acc: seq[PositionData] = @[]
      var storage = newSlidingDeque(20, proc (d: PositionData) = acc.add(d))
      
      # test legal record
      storage.record(0, "A", 'A')
      storage.record(1, "G", 'G')
      storage.record(2, "A", 'T')
      doAssert storage.deq.len == 3
      
      doAssert storage.flushUpTo(2) == 2
      
      doAssert storage.beginning == 2
      doAssert storage.deq.len == 1
      doAssert acc.len == 2

      doAssert acc[0].referenceIndex == 0
      doAssert acc[0].referenceBase == 'A'

      doAssert acc[1].referenceBase == 'G'
      doAssert acc[1].referenceIndex == 1

    block: # test flush without overlap
      var acc: seq[PositionData] = @[]
      var storage = newSlidingDeque(20, proc (d: PositionData) = acc.add(d))
      storage.record(0, "A", 'A')
      storage.record(1, "G", 'G')
      storage.record(2, "A", 'T')

      doAssert storage.flushUpTo(4) == 3

      doAssert storage.beginning == 4
      doAssert storage.deq.len == 0
      doAssert acc.len == 3

      doAssert acc[0].referenceIndex == 0
      doAssert acc[0].referenceBase == 'A'

      doAssert acc[1].referenceBase == 'G'
      doAssert acc[1].referenceIndex == 1

      doAssert acc[2].referenceBase == 'T'
      doAssert acc[2].referenceIndex == 2

    block:
      var actual : seq[PositionData] = @[]
      var storage = newSlidingDeque(20, proc (d: PositionData): void = actual.add(d))
      
      assert storage.flushUpTo(0) == 0
      storage.record(0,"A", 'A')
      storage.record(1,"A", 'A')
      storage.record(2, "C", 'A')
      storage.record(3,"A", 'A')
      doAssert storage.deq.len == 4
      doAssert storage.beginning == 0

      assert storage.flushUpTo(0) == 0
      storage.record(0,"A", 'A')
      storage.record(1,"T", 'A')
      storage.record(2, "G", 'A')
      storage.record(3,"A", 'A')
      doAssert storage.deq.len == 4
      doAssert storage.beginning == 0

      assert storage.flushUpTo(0) == 0
      storage.record(0,"A", 'A')
      storage.record(1,"T", 'A')
      storage.record(2, "-AC", 'A')
      storage.record(3,"G", 'G')
      storage.record(4,"C", 'G')
      doAssert storage.deq.len == 5, $storage.deq.len
      doAssert storage.beginning == 0

      assert storage.flushUpTo(10) == 5
      storage.record(10, "A", 'A')
      doAssert storage.deq.len == 1
      doAssert storage.beginning == 10

      var expected = @[
        {"A": 3}.toTable, 
        {"A": 1, "T": 2}.toTable, 
        {"C": 1, "G": 1, "-AC": 1}.toTable,
        {"A": 2, "G": 1}.toTable,
        {"C": 1}.toTable
      ]

      for idx, pair in zip(actual, expected):
        echo pair[0][]
        echo pair[1]
        doAssert pair[0].events == pair[1], $pair[0].events & " " & $pair[1]






