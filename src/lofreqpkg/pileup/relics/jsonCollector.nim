import json
import storage/containers/positionData
import interfaces/iCollector

type JsonCollector*[TNext] = ref object
  next: TNext

proc newJsonCollector*[TNext](next: TNext): JsonCollector[TNext] =
  JsonCollector[TNext](next: next)

proc submit*[TData](self: JsonCollector, data: TData): void =
  self.next.submit(%data)

proc getICollector*(self: JsonCollector): ICollector =
  return (
        submit: proc(data: PositionData): void =
          self.submit(data)
        )
