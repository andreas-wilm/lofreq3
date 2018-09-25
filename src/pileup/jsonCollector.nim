import json
import storage/containers/positionData
import interfaces/iCollector

type JsonCollector* = ref object
  chromosome: string


proc newJsonCollector*(chromosome: string): JsonCollector =
  JsonCollector(chromosome: chromosome)

proc addChromosome(data: JsonNode, chromosome: string): JsonNode =
  data["chromosome"] = %chromosome
  return data

proc submit*(self: JsonCollector, data: PositionData): void =
  writeLine(stdout, (%data).addChromosome(self.chromosome))


proc getICollector*(self: JsonCollector): ICollector =
  return (
      submit: proc(data: PositionData): void =
        self.submit(data),
     )
