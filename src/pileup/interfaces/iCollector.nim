import ../storage/containers/positionData

type 
  ICollector* = tuple[
    submit: proc(data: PositionData): void {.closure.},
  ]
