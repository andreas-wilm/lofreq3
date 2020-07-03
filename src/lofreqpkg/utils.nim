## LoFreq: helper functions
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License


import math
#import strutils
import times

const REF_SYMBOL_AT_INDEL_FW* = '-'
const REF_SYMBOL_AT_INDEL_RV* = '_'
const DEFAULT_BLANK_SYMBOL* = '*'# missing base symbol


proc refSymbolAtIndel*(reverse: bool): char =# evil hack to support strand for "non-indels"
  if reverse:
    REF_SYMBOL_AT_INDEL_RV
  else:
    REF_SYMBOL_AT_INDEL_FW


proc decodeQual*(encqual: char, offset: uint = 33): uint =
  return uint(ord(encqual))-offset


proc encodeQual*(qual: uint, offset: uint = 33): char =
  assert qual+offset < 255
  result = chr(qual+offset)


## brief convert error probability to phred quality
proc prob2qual*(e: float): Natural =
  # FIXME handle 0.0 with caught exception?
  assert e>=0.0
  if e.classify == fcZero:
    return high(Natural)
  return Natural(round(-10.0 * log10(e)))


## brief convert phred quality to error probability
proc qual2prob*(q: Natural): float =
  return pow(10.0, float(-q)/10.0)


## brief generate date string
proc dateStr*(): string =
  var t = getTime().local()
  result = t.format("yyyyMMdd")


# poor man's testing
template testblock*(message: string, body: untyped) =
  block:
    echo "Testing " & message
    body


when isMainModule:

  testblock "decodeQual":
    doAssert decodeQual('#') == 2
    
  testblock "encodeQual":
    doAssert encodeQual(2) == '#'

  testblock "prob2qual":
    doAssert prob2qual(0.05) == 13
    doAssert prob2qual(0.01) == 20
    doAssert prob2qual(0.001) == 30

  testblock "qual2prob":
    doAssert abs(qual2prob(13) - 0.05) < 0.005
    doAssert abs(qual2prob(20) - 0.01) < 0.001
    doAssert abs(qual2prob(30) - 0.001) < 0.0001

  echo "OK: all tests passed"