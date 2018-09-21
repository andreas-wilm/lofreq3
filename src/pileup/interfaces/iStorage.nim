type 
  IStorage* = tuple[
    record: proc(position: int,value: string,refBase: char): void {.closure.},
    flushUpTo: proc(position: int): int {.closure.},
    flushAll: proc(): int {.closure.}
  ]
