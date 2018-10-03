import json
import qualityHistogram

type OperationData*[T] = object
  histogram: QualityHistogram[T]
  reverse: int
  forward: int


proc initOperationData*[T](): OperationData[T] {.inline.} =
  OperationData[T](histogram: initQualityHistogram[T](), reverse: 0, forward: 0)


proc add*[T](self: var OperationData, bases: T, quality: int,
             reverse: bool) {.inline.} =
  self.histogram.add(bases, quality)
  if reverse:
    self.reverse.inc
  else:
    self.forward.inc

# If I just make one generic method, it doesn't work so I had
# to 'pattern match'
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
