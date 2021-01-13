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
import sequtils

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
  ## for debugging in pileup but need to be removed before calling
  
  let events = toSeq(self.keys())
  for ev in events:# can't use pairs() because we need qHist to be writable
    var qHist = self[ev]
    # del not supported for CountTable but setting to 0 should have the same effect:
    # https://stackoverflow.com/questions/59160984/remove-key-from-counttable-in-nim
    # still need to track zombie entries with eventCounts, because len() doesn't change
    # (CountTable bug; see above)
    var eventCounts = 0
    for qual in qHist.keys():
      if qual == -1:# ignore everything filtered (marked as q=-1 in pileup)
          self[ev].del(qual)
      else:
        inc(eventCounts, self[ev][qual])
    # Now delete zombie entries, i.e those with only filtered events
    if eventCounts == 0:
      self.del(ev)


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
    # pair[0] is the event, e.g. G or * or -
    buff[$pair[0]] = %pair[1]
  result.fields = buff


proc `%`*(table: CountTable[int]): JsonNode {.inline.} =
  result = newJObject()

  var buff = initOrderedTable[string, JsonNode]()
  # we are only sorting so that output is stable and can be used for testing
  var srttable = deepCopy(table)
  srttable.sort()
  for pair in srttable.pairs:
    # pair[0] is the quality
    buff[$pair[0]] = %pair[1]
  result.fields = buff

