## LoFreq: vcf routines
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License

import utils
import strformat

type Dp4* = object
  refForward*: int
  refReverse*: int
  altForward*: int
  altReverse*: int

type InfoField* = object
  af*: float
  sb*: int
  dp4*: Dp4

type Variant* = ref object
  chrom*: string
  pos*: int
  id*: string
  refBase*: char
  alt*: string
  qual*: int
  filter*: string
  info*: InfoField


proc `$`*(v: Variant): string =
  let dp4 = fmt"{v.info.dp4.refForward},{v.info.dp4.refReverse},{v.info.dp4.altForward},{v.info.dp4.altReverse}"
  var info = fmt"AF={v.info.af:.6f};SB={v.info.sb};DP4={dp4}"
  fmt("{v.chrom}\t{v.pos}\t{v.id}\t{v.refBase}\t{v.alt}\t{v.qual}\t{v.filter}\t{info}")# fmt() to get tabs


## brief create vcf header
# FIXME check conformity. refFa and src likely missing since only known in plp
proc vcfHeader*(src: string = "", refFa: string = ""): string =
  result = "##fileformat=VCFv4.2\n"
  result = result & "##fileDate=" & dateStr() & "\n"
  if len(src)>0:
    result = result & "##source=" & src & "\n"
  if len(refFa)>0:
    result = result & "##reference=" & refFa & "\n"

  result = result & """##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""


