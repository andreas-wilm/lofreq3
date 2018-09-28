## The module provides an implementation of a `PositionData` and its operations.
## 'PositionData' serves as a main data structure used in performing the pileup.
## It keeps all the information for a single position on the reference
import qualityHistogram 
import json

type
  PositionData* = ref object
    referenceIndex*: int
    referenceBase*: char
    matches: QualityHistogram[char]
    deletions: QualityHistogram[string]
    insertions: QualityHistogram[string]
    # chromosome - injected after pileup is done in order to save space and time


proc newPositionData*(referenceIndex: int, referenceBase: char) : PositionData =
  ## Constructs a new PositionData object keeping the data for the given position
  ## on the reference. The second argument should provide the base appearing on the said
  ## position.
  PositionData(
    referenceIndex: referenceIndex,
    referenceBase: referenceBase,
    matches: initQualityHistogram[char](),
    insertions: initQualityHistogram[string](),
    deletions: initQualityHistogram[string]()
  )


proc addMatch*(self: var PositionData, base: char, quality: int) =
  ## Accounts for a match on the position represented by this object.
  ## A match can either be a true match (same base as the reference)
  ## or a mismatch (base different from the reference) as long as the 
  ## base is present. A match is defined by its base and its quality.
  self.matches.add(base, quality)


proc addInsertion*(self: var PositionData, bases: string, quality: int) =
  ## Accounts for an insertion on the position represented by this object.
  ## An insertion consists of one or more bases not present on the reference.
  ## It is defined by its value (one or more bases) and its quality.
  self.insertions.add(bases, quality)


proc addDeletion*(self: var PositionData, bases: string, quality: int) =
  ## Accounts for a deletion the position represented by this object.
  ## A Deletion  consists of one or more missing bases (wrt. the reference).
  ## It is defined by its value (one or more bases) and its quality.
  self.deletions.add(bases, quality)


proc `%`*(self: PositionData): JsonNode =
  result = %{
    "referenceIndex": %self.referenceIndex,
    "referenceBase": %($self.referenceBase),
    "matches": %self.matches,
    "insertions": %self.insertions,
    "deletions": %self.deletions
  }
