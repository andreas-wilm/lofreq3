## The module contains the procedures used in processing the data after all
## information is collected.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

# standard
import json
# third party
# project specific
import storage/containers/positionData


proc toJson*(data: PositionData): JsonNode =
  ## Converts the given PositionData object into a JsonNode.
  %data


proc toJsonAndPrint*(data: PositionData): void =
  ## Converts the given PositionData object into a JsonNode.
  writeLine(stdout, %data)


proc doNothing*(data: PositionData): void =
  discard true

