import eventData
import json

type
  PositionData* = ref object
    referenceIndex*: int
    referenceBase*: char
    events*: EventData
    # chromosome: string

proc newPositionData* : PositionData =
  PositionData(referenceIndex: 0, referenceBase: '/', events: initEventData())

proc newPositionData*(referenceIndex: int, referenceBase: char) : PositionData =
  PositionData(
    referenceIndex: referenceIndex,
    referenceBase: referenceBase,
    events: initEventData()
  )

proc increment*(self: var PositionData, value: string) =
  self.events.increment(value)

proc `$`*(self: var PositionData): string =
  "(" & "referenceIndex: " & $self.referenceIndex & ", referenceBase: " & $self.referenceBase & ", events: " & $self.events & ")"
    
proc `%`*(self: PositionData): JsonNode =
  result = %{
    "referenceIndex": %self.referenceIndex,
    "referenceBase": %($self.referenceBase),
    "events": %self.events
  }