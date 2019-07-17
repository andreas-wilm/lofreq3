## The module contains the procedures used in processing the data after all
## information is collected.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

import storage/containers/positionData
import json


proc getJsonPropertyInjector*[TValue](key: string,
                                      value: TValue
                                      ): (proc (data:JsonNode): void) =
  ## Returns a procedure which takes a JsonNode and injects the given chromsome.
  ## name into it.
  return proc (data: JsonNode): void =
    data[key] = %value


proc print*[TData](data: TData) : void  =
  ## Prints the given data to the standard output.
  writeLine(stdout, $data)


proc toJson*(data: PositionData): JsonNode =
  ## Converts the given PositionData object into a JsonNode.
  %data

proc doNothing*(data: PositionData): void =
  discard true

