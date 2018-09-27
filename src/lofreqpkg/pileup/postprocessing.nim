## All the functions used in processing the data after all information is collected.

import storage/containers/positionData
import json

proc chromosomeInjector*(chromosome: string): (proc (data:JsonNode): JsonNode) =
  ## Returns a function which takes a JsonNode and injects the given chromsome
  ## name into it.
  return proc (data: JsonNode): JsonNode =
    data["chromosome"] = %chromosome
    return data

proc printOutput*[TData](data: TData) : void  =
  ## Prints the given data to the standard output
  writeLine(stdout, $data)

proc toJson*(data: PositionData): JsonNode =
  ## Converts the given PositionData object into a JsonNode
  %data

proc somethingStupid*(data: JsonNode): JsonNode =
  ## Does something stupid
  data["something"] = %"stupid" 
  return data

proc deleteDeletions*(data: JsonNode): JsonNode =
  ##  the deletion field from the provided JsonNode
  data.delete("deletions")
  return data




