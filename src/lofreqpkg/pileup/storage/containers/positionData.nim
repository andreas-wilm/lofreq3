## The module provides an implementation of a `PositionData` object and its
## procedure.  'PositionData' serves as a main data structure used in
## performing the pileup. It keeps all the information concerning a single
## position on the reference. This module does not have any authority over how
## it is used and which values does it store. It does not check the integrety
## or the semantics of the given data. While this certainly provides more
## flexibility, it also opens up the door for abuse. The user is expected to
## pass valid data. The procedures' docs provide guidelines for their usage.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License


import operationData
import json


type PositionData* = ref object
    ## The 'PositionData' object keeping all information concerning one parti-
    ## cular position on the reference.
    refIndex*: int
    refBase: char
    chromosome: string
    matches: OperationData[char]
    deletions: OperationData[string]
    insertions: OperationData[string]


proc newPositionData*(refIndex: int, refBase: char,
                      chromosome: string) : PositionData {.inline.} =
  ## Constructs a new PositionData object keeping the data for
  ## the given position on the reference.
  ## The position (wrt. the reference) is provided by the first argument.
  ## The second argument should provide the base appearing on the said
  ## position.
  ## The third argument must specify the name of the chromosome the
  ## reference sequence belongs to.
  PositionData(
    refIndex: refIndex,
    refBase: refBase,
    chromosome: chromosome,
    matches: initOperationData[char](),
    insertions: initOperationData[string](),
    deletions: initOperationData[string]()
  )


proc addMatch*(self: var PositionData, base: char, quality: int,
               reverse: bool) {.inline.} =
  ## Accounts for a match on the position represented by this object.  In the
  ## most common basic biological use case, a match should either be a true
  ## match (same base as the reference) or a mismatch (base different from the
  ## reference) as long as the base is present.  A match is defined by its
  ## base, its quality and its strand. Additionaly, missing positions (after a
  ## deletion) can also be saved as matches marked with a special character.
  ## This is something decided outside of this module. The procedure only
  ## ensures the value is accounted for and does not check the data's
  ## validity.
  self.matches.add(base, quality, reverse)


proc addInsertion*(self: var PositionData, bases: string, quality: int,
                   reverse: bool) {.inline.} =
  ## Accounts for an insertion on the position represented by this object.
  ## In the most common biological use case, an insertion consists of one or 
  ## more bases not present on the reference. It is defined by its value (one
  ## or more bases), its quality and its strand. The procedure only ensures the
  ## value is accounted for and does not check the data's validity.
  self.insertions.add(bases, quality, reverse)


proc addDeletion*(self: var PositionData, bases: string, quality: int,
                  reverse: bool) {.inline.} =
  ## Accounts for a deletion the position represented by this object.  In the
  ## most common biological use case, a deletion  consists of one or more
  ## missing bases (wrt. the reference). It is defined by its value (one or
  ## more bases), its quality and its strand. The procedure only ensures the
  ## value is accounted for and does not check the data's validity.
  self.deletions.add(bases, quality, reverse)


proc `%`*(self: PositionData): JsonNode {.inline.} =
  %{
    "chromosome": %self.chromosome,
    "refIndex": %self.refIndex,
    "refBase": %($self.refBase),
    "matches": %self.matches,
    "insertions": %self.insertions,
    "deletions": %self.deletions
  }
