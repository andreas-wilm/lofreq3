import sets
import sequtils
import tables
import containers/positionData 
import interfaces/iStorage 
import math

type SlidingTable* = ref object
  table: Table[int, PositionData]
  indices: HashSet[int]
  initialSize: int
  submit : proc (data: PositionData): void

proc newSlidingTable*(initialSize: int, 
                     submitProc: proc (data: PositionData): void
                    ): SlidingTable =
  let adjustedSize = nextPowerOfTwo(initialSize)
  SlidingTable(
    table: initTable[int, PositionData](adjustedSize),
    indices: initSet[int](),
    initialSize: adjustedSize, 
    submit: submitProc
  )

proc submitTable(self: SlidingTable, table: var Table[int, PositionData]): void =
  # todo implement
  # asyncronous function
  # Submits the current deque to another thread for handling
  for value in table.values:
    self.submit(value)

proc resetTable(self: SlidingTable) =
  self.submitTable(self.table)
  self.table = initTable[int, PositionData](self.initialSize)

proc flushUpTo*(self: SlidingTable, current: int): int =
  # try to turn this into a regular filter
  iterator filterSet[T](s: HashSet, predicate: proc (x: T): bool): T = 
    for element in s:
      if predicate(element):
        yield element

  for index in self.indices.filterSet(proc (i: int): bool = i < current):
    var value : PositionData
    discard self.table.take(index, value)
    self.submit(value) # make async
    result.inc

proc flushAll*(self: SlidingTable): int =
  result = self.table.len
  self.resetTable()
  self.indices = initSet[int]()

proc record*(self: SlidingTable, position: int, value: string, 
                   refBase: char): void =
  self.indices.incl(position)
  self.table.mgetOrPut(position, newPositionData(position, refBase)).increment(value)

proc recordStart(self: SlidingTable, position: int, value: string,
                 refBase: char): void =
  discard self.flushUpTo(position) # make async
  self.record(position, value, refBase)

proc getIStorage*(self: SlidingTable): IStorage =
  return (
        record: proc(position: int, value: string, refBase: char): void = 
          self.record(position, value, refBase),
        flushUpTo: proc(position: int): int = 
          self.flushUpTo(position),
        flushAll: proc(): int  = 
          self.flushAll()
       )


