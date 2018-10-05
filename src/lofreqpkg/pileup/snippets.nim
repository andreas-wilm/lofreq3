proc chromosomeInjector*(chromosome: string): (proc (data:JsonNode): JsonNode) =
  return proc (data: JsonNode): JsonNode =
    data["chromosome"] = %chromosome
    return data

proc toJson(data: PositionData): JsonNode =
  %data

proc addChromosome(data: JsonNode): JsonNode =
  data["chromosome"] = %"ref"
  return data

proc removeDeletions(data: JsonNode): JsonNode =
  data.delete("deletions");
  return data

proc somethingStupid(data: JsonNode): JsonNode =
  data["something"] = %"stupid" 
  return data

proc addChromosomeVoid(data: JsonNode): void =
  data["chromosome"] = %"ref"

proc removeDeletionsVoid(data: JsonNode): void =
  data.delete("deletions");

proc somethingStupidVoid(data: JsonNode): void =
  data["something"] = %"stupid" 

proc printIt[T](data: T): void =
  echo data

proc printStringLength[T](data: T): void =
  echo ($data).len


  var handler = (proc (x: PositionData): JsonNode = %x)
    .thenDo(proc (x:JsonNode): void = x["chromosome"] = %name)
    .then(proc (x:JsonNode): void = echo x)

  var handler = toJson
    .then(addChromosome)
    .then(somethingStupid)
    .thenDo(printStringLength)
    .then(removeDeletions)
    .thenDo(printIt)
    .done()

  var storage = newSlidingDeque(200, handler)
  

