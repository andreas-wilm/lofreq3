import hts

type 
  ISequence* =  tuple[
    baseAt: proc(index: int): char {.closure.},
    substring: proc(first, last: int): string {.closure.},
  ]

proc getISequence*(fai: Fai): ISequence =
  let sequence = fai.get(fai[0])
  return (
      baseAt: proc (index: int): char = sequence[index],
      substring: proc (first, last: int): string = sequence[first..last]
    )