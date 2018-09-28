import json

type ChromosomeInjector*[TNext] = ref object
  chromosome: string
  next: TNext

proc newChromosomeInjector*[TNext](next: TNext,
                                   chromosome: string): ChromosomeInjector[TNext] =
  ChromosomeInjector[TNext](chromosome: chromosome, next: next)

proc submit*(self: ChromosomeInjector, data: JsonNode) : void =
  data["chromosome"] = %self.chromosome
  self.next.submit(data)
