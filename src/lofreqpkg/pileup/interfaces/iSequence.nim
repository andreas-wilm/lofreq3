import hts

type 
  ISequence* =  tuple[
    baseAt: proc(index: int): char {.closure.},
    substring: proc(first, last: int): string {.closure.},
    len: int
  ]


proc loadSequence*(fai: Fai, name: string): ISequence =
  if fai.isNil:
    # fake all N sequence (length -1!)
    return (
      baseAt: proc (index: int): char {.closure.} = 'N',
      substring: proc (first, last: int): string {.closure.} = "N",
      len: -1
    )
  else:
    let sequence = fai.get(name)
    return (
      baseAt: proc (index: int): char = sequence[index],
      substring: proc (first, last: int): string = sequence[first..last],
      len: sequence.len
    )    
