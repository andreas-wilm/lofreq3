import json
import storage/containers/positionData

proc createJsonMessage*(data: PositionData): string = $(%data) & "\c\l"