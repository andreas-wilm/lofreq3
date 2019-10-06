## LoFreq: helper functions
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License


import math
import strutils
import times

proc decodeASCIIQual*(qual: char,  offset: Natural = 33): Natural =
  return ord(qual)-offset


## brief convert error probability to phred quality
proc prob2qual*(e: float): Natural =
  # FIXME handle 0.0 with caught exception?
  assert e>=0.0
  if e.classify == fcZero:
    return high(Natural)
  return Natural(-10.0 * log10(e))


## brief convert phred quality to error probability
proc qual2prob*(q: Natural): float =
  return pow(10.0, float(-q)/10.0)


## brief generate date string
proc dateStr*(): string =
  var t = getTime().local()
  result = t.format("yyyyMMdd")

