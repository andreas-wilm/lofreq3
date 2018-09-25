import eventData
import json

type
  PositionData* = ref object
    referenceIndex*: int
    referenceBase*: char
    matches: EventData[char]
    deletions: EventData[string]
    insertions: EventData[string]
    # chromosome - injected after pileup is done in order to save space and time


proc newPositionData*(referenceIndex: int, referenceBase: char) : PositionData =
  PositionData(
    referenceIndex: referenceIndex,
    referenceBase: referenceBase,
    matches: initEventData[char](),
    insertions: initEventData[string](),
    deletions: initEventData[string]()
  )


proc addMatch*(self: var PositionData, base: char, quality: int) =
  self.matches.add(base, quality)


proc addInsertion*(self: var PositionData, bases: string, quality: int) =
  self.insertions.add(bases, quality)


proc addDeletion*(self: var PositionData, bases: string, quality: int) =
  self.deletions.add(bases, quality)


proc `%`*(self: PositionData): JsonNode =
  result = %{
    "referenceIndex": %self.referenceIndex,
    "referenceBase": %($self.referenceBase),
    "matches": %self.matches,
    "insertions": %self.insertions,
    "deletions": %self.deletions
  }