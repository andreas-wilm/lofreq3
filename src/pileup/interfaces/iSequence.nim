import hts

type 
  ISequence* =  tuple[
    baseAt: proc(index: int): char {.closure.},
    substring: proc(first, last: int): string {.closure.},
  ]

proc getISequence*(fai: Fai): ISequence = 
  let name = fai[0]
  return (
      baseAt: proc (index: int): char = fai.get(name, index, index)[0],
      substring: proc (first, last: int): string = fai.get(name, first, last)
    )