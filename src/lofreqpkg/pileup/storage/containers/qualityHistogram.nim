## The module provides an implementation of the 'QualityHistogram' object and
## its methods. 'QualityHistogram' should be used to count all events on the
## same position and of the same kind. Events are uniquely determined by their
## position on the reference, their operation kind, their value, their quality
## and their strand. This object and its procedures distinguish entries based
## only on their values and qualities. Because different kinds of operations
## have different types for their values, the object and its methods are
## parameterized.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

import tables
import json


## Defines a 'QualityHistogram' type with a type parameter specifying the type
## of the event value.
type QualityHistogram*[T] = Table[T, CountTable[int]]


func initQualityHistogram*[T](): QualityHistogram[T] {.inline.} =
  ## Constructs a new QualityHistogram object. The type parameter T sets the
  ## type of event values.
  initTable[T, CountTable[int]]()


proc coverage*[T](table: QualityHistogram[T]): Natural =
  var s = 0
  for p1 in table.pairs:
    #s = p1[1]
    for p2 in p1[1].pairs:
      s += p2[1]
  return s


proc clean*[T](self: var QualityHistogram[T]): void =
  ## removes filtered entries, i.e those with q<0 that are kept
  ## for debugging in pileup but not to be removed before calling

  for event in self.keys():# can't use pairs because we need qHist to be writable
    var qHist = self[event]
    # del not supported for CountTable but setting to 0 should have the same effect:
    # https://stackoverflow.com/questions/59160984/remove-key-from-counttable-in-nim
    # still need to track zombie entries with eventCounts, because len() doesn't change
    # (CountTable bug; see above)
    var eventCounts = 0
    for qual in qHist.keys():
      if qual == -1:# ignore everything filtered (marked as q=-1 in pileup)
        self[event][qual] = 0# qHist is just a copy here!
        # in the future we can use this:self[event].del(qual)
      else:
        inc(eventCounts, self[event][qual])
    # delete zombie entries (those with only filtered events)
    #if len(self[event]) == 0: once this is fixed in Nim we can remove eventCounts
    if eventCounts == 0:
      self.del(event)


proc set*[T](self: var QualityHistogram[T], value: T,
             quality: int, count: int): void {.inline.} =
  ## Accounts for a event with the given value and the given quality.
  discard self.hasKeyOrPut(value, initCountTable[int]())
  self[value][quality] = count


proc add*[T](self: var QualityHistogram[T], value: T,
             quality: int): void {.inline.} =
  ## Accounts for a event with the given value and the given quality.
  discard self.hasKeyOrPut(value, initCountTable[int]())
  self[value].inc(quality)


proc `%`*[T](table: var QualityHistogram[T]): JsonNode {.inline.} =
  result = newJObject()

  var buff = initOrderedTable[string, JsonNode]()
  for pair in table.pairs:
    buff[$pair[0]] = %pair[1]

  result.fields = buff


proc `%`*(table: CountTable[int]): JsonNode {.inline.} =
  result = newJObject()

  var buff = initOrderedTable[string, JsonNode]()
  for pair in table.pairs:
    buff[$pair[0]] = %pair[1]

  result.fields = buff

