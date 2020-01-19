## The module implements a genomics region
##
## - Author: Andreas Wilm <andreas.wilm@gmail.com>
## - License: The MIT License

import nre
import strformat
import strutils

# FIXME should really be uints or int64
# zero based, half open as bed
type Region* = object
  sq*: string
  s*: Natural
  e*: Natural


#proc autoFillRegion*(reg: Region, targets: seq[Target]): Region =
#  result = reg
#  var foundTarget = false
#  for t in targets:
#    if t.name == reg.sq:
#      result.s = 0
#      result.e = (int)t.length
#      foundTarget = true
#      break
#  if not foundTarget:
#    raise newException(ValueError,
#      fmt"Couldn't find {reg.sq} in BAM header. Valid entries are {targets}")


proc regionFromStr*(regStr: string): Region =
  ## expects region string in the form of seq[:start-end],
  ## where start-end are inclusive and 1-based (and >0).
  ## indices are stored as zero based, half open interval.
  ## if start and end  are not given, both will be set to 0
  var matches = regStr.match(re"(.+):(\d+)-(\d+)")
  if matches.isSome:
    # chr:start-end
    result.sq = matches.get.captures[0]
    # FIXME test correct boundaries
    result.s = matches.get.captures[1].parseInt-1
    result.e = matches.get.captures[2].parseInt
    doAssert result.s >= 0 and result.s < result.e# FIXME rather raise exception
  else:
    # chrom only
    doAssert regStr.contains(':') == false# FIXME proper error message
    result.sq = regStr
    # s=e means no indices given
    result.s = 0
    result.e = 0
  doAssert len(result.sq)>0


iterator bedParser(bedFile: string): Region =
  var reg: Region
  for line in lines(bedFile):
    if line.startswith("#") or  line.startswith("track"):
      continue
    var fields = line.strip().split('\t', 5)
    if len(fields) < 3:
      raise
    reg.sq = fields[0]
    reg.s = parseInt(fields[1])
    reg.e = parseInt(fields[2])
    doAssert reg.e > reg.s and reg.s >= 0# FIXME duplication
    yield reg


iterator parseRegionsStr(regionsStr: string): Region =
  for regStr in regionsStr.split(','):
    var reg = regionFromStr(regStr)
    doAssert reg.e > reg.s and reg.s >= 0# FIXME duplication
    #if reg.s == 0 and reg.e == 0:# only sq given instead of full region
    #var targets = targets(bam.hdr)
    #reg = auto_fill_region(reg, targets)
    yield reg


iterator getRegions*(regionsStr = "", bedFile = ""): Region =
  # allows in theory to use both regionsStr and bedFile
  if len(regionsStr)>0:
    for reg in parseRegionsStr(regionsStr):
      yield reg
  if len(bedFile)>0:
     for reg in bedParser(bedFile):
       yield reg


proc `$`*(r: Region): string =
  fmt"{r.sq}:{r.s+1}-{r.e}"
