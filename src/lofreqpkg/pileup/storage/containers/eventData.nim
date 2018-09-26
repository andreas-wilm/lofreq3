## The module provides an implementation of the 'EventData' object and its methods.
## 'EventData' should be used to count all events on the same position and of the same
## kind. Events are uniquely determined by their position on the reference,
## their kind, their value and their quality. Because different kinds of events have
## different types for their values, the object and its methods are parameterized.

import tables
import json

type EventData*[T] = Table[T, CountTable[int]]


func initEventData*[T](): EventData[T] =
  ## Constructs a new EventData object. The type parameter T sets the
  ## type of event values.
  initTable[T, CountTable[int]]()


proc add*[T](self: var EventData[T], value: T, quality: int): void =
  ## Accounts for a event with a given value and quality.
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

