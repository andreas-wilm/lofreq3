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
import tables
import strutils


type OperationData*[T] = object
  ## The 'OperationData' type. It makes and provides a histogram on operation
  ## values and their qualities
  histogram: QualityHistogram[T]


proc initOperationData*[T](): OperationData[T] {.inline.} =
  ## Creates a new 'OperationData' object. All that it needs is the type for
  ## the operation values.
  OperationData[T](histogram: initQualityHistogram[T]())


proc add*[T](self: var OperationData, bases: T, quality: int,
             reverse: bool): void {.inline.} =
  ## Accounts for one specific operation. Distinct oprations of the same kind
  ## are determined by their value, their quality and their strand. The
  ## histogram, however, does not take the strand into the account. In this
  ## case, distinct operations are determined only by their base and their
  ## quality. The strand information is counted separately.

  if reverse:
    self.histogram.add(bases.toLowerAscii(), quality)
  else:
    self.histogram.add(bases, quality)


# If I just make one generic method, it doesn't work so I had to 'pattern
# match'.
proc `%`*(self: var OperationData[char]): JsonNode {.inline.} =
  %self.histogram


proc `%`*(self: var OperationData[string]): JsonNode {.inline.} =
  %self.histogram
