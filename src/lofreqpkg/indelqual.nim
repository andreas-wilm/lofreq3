## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <andreas.wilm@gmail.com>
## - License: The MIT License

# standard
import tables
import strutils
import strformat

# third party
import hts

# project specific
import utils

const DINDELQ = "!MMMLKEC@=<;:988776"# 1-based 18
# ? const DINDELQ2 = "!CCCBA;963210/----,"#  *10 

const BI_TAG = "BI"
const BD_TAG = "BD"


proc skipRead(rec: Record): bool =
  if rec.flag.secondary or rec.flag.qcfail or rec.flag.dup or rec.flag.unmapped:
    return true
  return false


 # a groupby variant that splits a string into runs of identical
 # characters which is yielded included the length of the run
iterator groupBy(query: string): (int, char) =
  doAssert len(query) != 0
  var runLength = 1 
  var runChar = query[0]
  for c in query[1 .. ^1]:
    if c == runChar:
      inc runLength
    else:
      yield (runLength, runChar)
      runLength = 1 
      runChar = c
  yield (runLength, runChar)


proc findHomopolymerRuns(query: string): seq[int] =
  for (l, c) in groupBy(query):
    result.add(l)
    for i in countup(1, l-1):
      result.add(1)
      

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


proc getDindelQual(rec: Record, homopolymerRuns: seq[int]): (string, string) = 
  var rpos = rec.start;# coordinate on reference x
  let rlen = len(homopolymerRuns)
  var dindelq: string
  for ce in rec.cigar:
    let co = ce.op
    let cl = ce.len
    if query(consumes(ce)) and reference(consumes(ce)):#M=X
      for i in 0 ..< cl:
        if rpos > rlen-2:
          dindelq.add(DINDELQ[0])
        else:
          if homopolymerRuns[rpos+1] > 18:
            dindelq.add(DINDELQ[0])
          else:
            dindelq.add(DINDELQ[homopolymerRuns[rpos+1]])
          inc rpos
    elif reference(consumes(ce)):#DN
      rpos += cl
    elif query(consumes(ce)):#IS
      dindelq.add(DINDELQ[0].repeat(cl))

  when not defined(release):
    var query: string
    discard rec.sequence(query)
    if len(dindelq) != len(query):
      stderr.writeLine("DEBUG dindelq=" & dindelq)
      stderr.writeLine("DEBUG   query=" & query)
      stderr.writeLine("DEBUG   cigar=" & $rec.cigar)
  assert len(dindelq) == len(query)
      
  result = (dindelq, dindelq)


proc parseIndelArg(arg: string): (char, char) =
  var args = split(arg, ",")
  let iq = parseUInt(args[0])
  var dq: uint
  if len(args) == 1:
    dq = iq
  elif len(args) == 2:
    dq = parseUInt(args[1])
  else:
    raise newException(ValueError,
      fmt"Couldn't find parse indel quality arg {arg}. Should be either indelq or insq,delq")
  result = (encodeQual(iq), encodeQual(dq))


proc indelqual*(faFname: string, bamInFname: string, uniform: string = "") =

  var fai: Fai
  # keeping all observed reference homopolymers in memory
  # this should only be needed if file isn't sorted. FIXME warn?
  var homopolymersPerChrom = initTable[string, seq[int]]()
  var iBam: Bam
  #var oBam: Bam
  var insQual: char
  var delQual: char
  if len(uniform) > 0:
    (insQual, delQual) = parseIndelArg(uniform)
 
  if not open(fai, faFname):
    quit("Could not open reference sequence file: " & faFname)

  open(iBam, bamInFname, fai=faFname)
  
  stdout.write($iBam.hdr)
  #assert outFormat in @["BAM", "SAM", "CRAM"]
  #open(oBam, "-", mode="w" & outFormat, fai=faFname)
  #obam.write_header(iBam.hdr)

  var bi: string
  var bd: string
  for rec in iBam:
    var chrom = rec.chrom

    if skipRead(rec):
      #obam.write(rec)
      echo $rec.tostring()
      continue

    # for dindel: load reference is needed, compute homopolymers and set bi and bd 
    if len(uniform) == 0:
      if not homopolymersPerChrom.hasKey(chrom):
        let refsq = fai.get(chrom)
        homopolymersPerChrom[chrom] = findHomopolymerRuns(refsq)
      (bi, bd) = getDindelQual(rec, homopolymersPerChrom[chrom])
        
    else:
      let l = rec.b.core.l_qseq# no len function in htsnim?
      bi = insQual.repeat(l)
      bd = delQual.repeat(l)

    echo createRec(rec, bi, bd)
    #obam.write(createRec(rec, bi, bd))

    #oBam.close()


when isMainModule:
  testblock "findHomopolymerRuns":
    let x = findHomopolymerRuns("AACCCTTTTA")
    doAssert x == @[2, 1, 3, 1, 1, 4, 1, 1, 1, 1]

  testblock "parseIndelArg":
    doAssert ('#', '#') == parseIndelArg("2")
    doAssert ('#', '#') == parseIndelArg("2,2")
    doAssert ('5', '6') == parseIndelArg("20,21")


  echo "OK: all tests passed"