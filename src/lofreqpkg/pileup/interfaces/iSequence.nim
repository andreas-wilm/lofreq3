import hts

type 
  ISequence* =  tuple[
    baseAt: proc(index: int): char {.closure.},
    substring: proc(first, last: int): string {.closure.},
    len: int
  ]

proc loadSequence*(fai: Fai, name: string): ISequence =
  let sequence = fai.get(name)
  return (
      baseAt: proc (index: int): char = sequence[index],
      substring: proc (first, last: int): string = sequence[first..last],
      len: sequence.len
    )
