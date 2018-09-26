import json
type StreamWriter* = ref object

proc newStreamWriter*(): StreamWriter =
  StreamWriter()

proc submit3*[TData](self: StreamWriter, data: TData) : void =
  writeLine(stdout, $data)
