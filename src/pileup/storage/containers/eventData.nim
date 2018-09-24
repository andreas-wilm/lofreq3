import tables
import json

type EventData*[T] = Table[T, CountTable[int]]


func initEventData*[T](): EventData[T] = 
  initTable[T, CountTable[int]]()


proc add*[T](self: var EventData[T], value: T, quality: int): void =
  discard self.hasKeyOrPut(value, initCountTable[int]())
  self[value].inc(quality)


proc `%`*[T](table: var EventData[T]): JsonNode =
  result = newJObject()
  
  var buff = initOrderedTable[string, JsonNode]()
  for pair in table.pairs:
    buff[$pair[0]] = %pair[1]

  result.fields = buff


proc `%`*(table: CountTable[int]): JsonNode =
  result = newJObject()
  
  var buff = initOrderedTable[string, JsonNode]()
  for pair in table.pairs:  
    buff[$pair[0]] = %pair[1]

  result.fields = buff

