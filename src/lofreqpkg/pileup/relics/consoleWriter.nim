import json
type StreamWriter* = ref object

proc newStreamWriter*(): StreamWriter =
  StreamWriter()

proc submit*[TData](self: StreamWriter, data: TData) : void =
  writeLine(stdout, $data)
