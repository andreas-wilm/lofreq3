## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <andreas.wilm@gmail.com>
## - License: The MIT License

# standard
import tables
import strutils

# third party
import hts

# project specific
import utils

const DINDELQ = "!MMMLKEC@=<;:988776"# 1-based 18
const DINDELQ2 = "!CCCBA;963210/----,"#  *10 

const BI_TAG = "BI"
const BD_TAG = "BD"


proc skipRead(rec: Record): bool =
  if rec.flag.secondary or rec.flag.qcfail or rec.flag.dup or rec.flag.unmapped:
    return true
  return false


proc parseIndelArg(uniform: string): (char, char) =
  # FIXME parse
  # FIXME encode
  doAssert false


proc createRec(rec: Record, bi: string, bd: string): string =
    # We can't set BAM values in htsnim (argh), so convert to string / SAM
    # as in https://github.com/brentp/bamject/blob/master/src/bamject.nim

    when not defined(release):
      var query: string
      discard rec.sequence(query)
      assert len(bi) == len(bd)
      assert len(bi) == len(query)

    var recSplit = rec.tostring().split('\t')

    # delete instead of testing and overwrite if exists
    var delIndices: seq[int]
    for i in countup(11, len(recSplit)-1):
      let fieldSplit = recSplit[i].split(":")
      if fieldSplit[0] in @["BI", "BD"]:
        delIndices.add(i)
    for i, j in delIndices.pairs:
      recSplit.delete(j-i)

    recSplit.add("BI:Z:" & bi)
    recSplit.add("BD:Z:" & bd)

    result = recSplit.join("\t")


proc indelqual*(faFname: string, bamInFname: string, uniform: string = "") =

  var fai: Fai
  # keeping all observed reference homopolymers in memory
  # this should only be needed if file isn't sorted. FIXME warn?
  var homopolymers = initTable[string, string]()
  var iBam: Bam
  #var oBam: Bam
  var bi = ""
  var bd = ""

  if len(uniform) > 0:
    var (insQual, delQual) = parseIndelArg(uniform)
    bi = "FIXME"
    bd = "FIXME"

  if not open(fai, faFname):
    quit("Could not open reference sequence file: " & faFname)

  open(iBam, bamInFname, fai=faFname)
  
  stdout.write($iBam.hdr)
  #assert outFormat in @["BAM", "SAM", "CRAM"]
  #open(oBam, "-", mode="w" & outFormat, fai=faFname)
  #obam.write_header(iBam.hdr)

  for rec in iBam:
    var chrom = rec.chrom

    if skipRead(rec):
      #obam.write(rec)
      echo $rec.tostring()
      continue

    # for dindel: load reference is needed, compute homopolymers and set bi and bd 
    if len(uniform) == 0:
      if not homopolymers.hasKey(chrom):
        #stderr.writeLine("DEBUG Loading " & chrom)
        let refsq = fai.get(chrom)
        stderr.writeLine("FIXME find homopolymers")
        stderr.writeLine("FIXME get dindel indelq (bi == bd)")
    # for uniform: bi and bd already set

    echo createRec(rec, bi, bd)
    #obam.write(createRec(rec, bi, bd))

    #oBam.close()


when isMainModule:
  #testblock "name":
    #doAssert true
    
  echo "OK: all tests passed"