## The module implements a data structure used to store all information
## relevant to a particular kind of an alignment operation. Traditionally, the
## possible operations are matches, insertions and deletions. 'OperationData'
## is a paremeterized type because different kinds of operations are defined in
## terms of different types (e.g. single 'char' for matches, 'string' for
## insertions and deletions). The module also provides a '%' (toJson) procedure
## which converts an 'OperationData' object into a JsonNode.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

import json
import qualityHistogram


type OperationData*[T] = object
  ## The 'OperationData' type. It makes and provides a histogram on operation
  ## values and their qualities, as well as a total number of forward and
  ## reverse operations.
  histogram: QualityHistogram[T]
  reverse: int
  forward: int


proc initOperationData*[T](): OperationData[T] {.inline.} =
  ## Creates a new 'OperationData' object. All that it needs is the type for
  ## the operation values.
  OperationData[T](histogram: initQualityHistogram[T](),
                   reverse: 0, forward: 0)


proc add*[T](self: var OperationData, bases: T, quality: int,
             reverse: bool): void {.inline.} =
  ## Accounts for one specific operation. Distinct oprations of the same kind
  ## are determined by their value, their quality and their strand. The
  ## histogram, however, does not take the strand into the account. In this
  ## case, distinct operations are determined only by their base and their
  ## quality. The strand information is counted separately.
  self.histogram.add(bases, quality)
  if reverse:
    self.reverse.inc
  else:
    self.forward.inc


# If I just make one generic method, it doesn't work so I had to 'pattern
# match'.
proc `%`*(self: var OperationData[char]): JsonNode {.inline.} =
  %{
    "histogram": %self.histogram,
    "reverse": %self.reverse,
    "forward": %self.forward
  }


proc `%`*(self: var OperationData[string]): JsonNode {.inline.} =
  %{
    "histogram": %self.histogram,
    "reverse": %self.reverse,
    "forward": %self.forward
  }
