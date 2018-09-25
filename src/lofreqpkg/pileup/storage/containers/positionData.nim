import eventData
import tables
import json


type
  PositionData* = ref object
    referenceIndex*: int
    referenceBase*: char
    matches: Table[char, CountTable[int]]
    deletions: Table[string, CountTable[int]]
    insertions: Table[string, CountTable[int]]
    # chromosome: string

proc newPositionData* : PositionData =
  PositionData(
    referenceIndex: 0, 
    matches: initTable[char, CountTable[int]](),
    insertions: initTable[string, CountTable[int]](),
    deletions: initTable[string, CountTable[int]]()
  )

proc newPositionData*(referenceIndex: int, referenceBase: char) : PositionData =
  PositionData(
    referenceIndex: referenceIndex,
    referenceBase: referenceBase,
    matches: initTable[char, CountTable[int]](),
    insertions: initTable[string, CountTable[int]](),
    deletions: initTable[string, CountTable[int]]()
  )

proc addMatch*(self: var PositionData, base: char, quality: int) =
  discard self.matches.hasKeyOrPut(base, initCountTable[int]())
  self.matches[base].inc(quality)


proc addInsertion*(self: var PositionData, bases: string, quality: int) =
  discard self.insertions.hasKeyOrPut(bases, initCountTable[int]())
  self.insertions[bases].inc(quality)

proc addDeletion*(self: var PositionData, bases: string, quality: int) =
  discard self.deletions.hasKeyOrPut(bases, initCountTable[int]())
  self.deletions[bases].inc(quality)


proc `$`*(self: var PositionData): string =
  "(" & "referenceIndex: " & $self.referenceIndex & ", referenceBase: " & $self.referenceBase & ", events: " & $self.matches & ")"

proc `%`*(table: CountTable[int]): JsonNode =
  result = newJObject()
  
  var buff = initOrderedTable[string, JsonNode]()
  for pair in table.pairs:  
    buff[$pair[0]] = %pair[1]

  result.fields = buff

proc `%`*[T](table: Table[T, CountTable[int]]): JsonNode =
  result = newJObject()
  
  var buff = initOrderedTable[string, JsonNode]()
  for pair in table.pairs:
    buff[$pair[0]] = %pair[1]

  result.fields = buff

proc `%`*(self: PositionData): JsonNode =
  result = %{
    "referenceIndex": %self.referenceIndex,
    "referenceBase": %($self.referenceBase),
    "matches": %self.matches,
    "insertions": %self.insertions,
    "deletions": %self.deletions
  }