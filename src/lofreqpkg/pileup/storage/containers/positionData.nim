## The module provides an implementation of a `PositionData` and its operations.
## 'PositionData' serves as a main data structure used in performing the pileup.
## It keeps all the information for a single position on the reference
import operationData
import json

type PositionData* = ref object
    referenceIndex*: int
    referenceBase*: char
    chromosome: string
    matches: OperationData[char]
    deletions: OperationData[string]
    insertions: OperationData[string]


proc newPositionData*(referenceIndex: int, referenceBase: char,
                      chromosome: string) : PositionData {.inline.} =
  ## Constructs a new PositionData object keeping the data for 
  ## the given position on the reference.
  ## The second argument should provide the base appearing on the said
  ## position.
  PositionData(
    referenceIndex: referenceIndex,
    referenceBase: referenceBase,
    chromosome: chromosome,
    matches: initOperationData[char](),
    insertions: initOperationData[string](),
    deletions: initOperationData[string]()
  )


proc addMatch*(self: var PositionData, base: char, quality: int,
               reverse: bool) {.inline.} =
  ## Accounts for a match on the position represented by this object.
  ## A match can either be a true match (same base as the reference)
  ## or a mismatch (base different from the reference) as long as the 
  ## base is present. A match is defined by its base and its quality.
  self.matches.add(base, quality, reverse)


proc addInsertion*(self: var PositionData, bases: string, quality: int,
                   reverse: bool) {.inline.} =
  ## Accounts for an insertion on the position represented by this object.
  ## An insertion consists of one or more bases not present on the reference.
  ## It is defined by its value (one or more bases) and its quality.
  self.insertions.add(bases, quality, reverse)


proc addDeletion*(self: var PositionData, bases: string, quality: int,
                  reverse: bool) {.inline.} =
  ## Accounts for a deletion the position represented by this object.
  ## A Deletion  consists of one or more missing bases (wrt. the reference).
  ## It is defined by its value (one or more bases) and its quality.
  self.deletions.add(bases, quality, reverse)


proc `%`*(self: PositionData): JsonNode {.inline.} =
  %{
    "chromosome": %self.chromosome,
    "referenceIndex": %self.referenceIndex,
    "referenceBase": %($self.referenceBase),
    "matches": %self.matches,
    "insertions": %self.insertions,
    "deletions": %self.deletions
  }
