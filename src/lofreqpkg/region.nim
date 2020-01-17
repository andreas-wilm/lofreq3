## The module implements a genomics region
##
## - Author: Andreas Wilm <andreas.wilm@gmail.com>
## - License: The MIT License

import nre
import strformat
import strutils

# FIXME should really be uints
type Region* = object
  sq*: string
  s*: int
  e*: int


proc `$`*(r: Region): string =
  fmt"{r.sq}:{r.s+1}-{r.e}"


proc reg_from_str*(regStr: string): Region =
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
  
