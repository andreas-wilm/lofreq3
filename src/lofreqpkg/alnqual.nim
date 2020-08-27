## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <andreas.wilm@gmail.com>
## - License: The MIT License

 
stderr.writeLine("WARNING: Hardcoded htslib path")

# standard
import tables
import strutils
import strformat

# third party
import hts

# project specific
import utils
import bam_md_ext

const AI_TAG = "ai"
const AD_TAG = "ad"
const BAQ_TAG = "lb"


proc skipRead(rec: Record): bool =
  if rec.flag.secondary or rec.flag.qcfail or rec.flag.dup or rec.flag.unmapped:
    return true
  return false


proc createRec(rec: Record, ai: string, ad: string): string =
    # We can't set BAM values in htsnim (argh), so convert to string / SAM
    # as in https://github.com/brentp/bamject/blob/master/src/bamject.nim

    when not defined(release):
      var query: string
      discard rec.sequence(query)
      assert len(ai) == len(ad)
      assert len(ai) == len(query)

    var recSplit = rec.tostring().split('\t')

    # delete existing tags 
    var delIndices: seq[int]
    for i in countup(11, len(recSplit)-1):
      let fieldSplit = recSplit[i].split(":")
      if fieldSplit[0] in @[AI_TAG, AD_TAG]:
        delIndices.add(i)
    for i, j in delIndices.pairs:
      recSplit.delete(j-i)

    # and add again
    recSplit.add(AI_TAG & ":Z:" & ai)
    recSplit.add(AD_TAG & ":Z:" & ad)

    result = recSplit.join("\t")


proc alnqual*(faFname: string, bamInFname: string) =

  var fai: Fai
  var iBam: Bam
  #var oBam: Bam
  # keeping all observed reference sequences in memory for speedup
  var refs = initTable[string, string]()

  if not open(fai, faFname):
    quit("Could not open reference sequence file: " & faFname)

  open(iBam, bamInFname, fai=faFname)
  
  stdout.write($iBam.hdr)
  #assert outFormat in @["BAM", "SAM", "CRAM"]
  #open(oBam, "-", mode="w" & outFormat, fai=faFname)
  #obam.write_header(iBam.hdr)

  # ../lofreq.git/src/lofreq/kprobaln_ext.c kpa_ext_glocal
  # ../lofreq.git/src/lofreq/bam_md_ext.c idaq and bam_prob_realn_core_ext

  for rec in iBam:
    #if has_ins or has_del:
    # FIXME then what?
    # delete existing tags
    # skip if not aligned

    # load reference if not cached
    var chrom = rec.chrom
    if not refs.hasKey(chrom):
      #stderr.writeLine("DEBUG Loading " & chrom)
      refs[chrom] = fai.get(chrom)

    #echo createRec(rec, bi, bd)
    #obam.write(createRec(rec, bi, bd))
    const baq_flag = 1
    const baq_extended = 1
    const idaq_flag = 1
    var query: string
    discard rec.sequence(query)
    var ai_str = newString(len(query))
    var ad_str = newString(len(query))
    var baq_str = newString(len(query))

    var bam_lf: bam_lf_t
    var rc = bam_prob_realn_core_ext(addr bam_lf, refs[chrom], 
                            baq_flag, baq_extended, idaq_flag, 
                            baq_str, ai_str, ad_str)
    notImplementedError

  #oBam.close()
  
when isMainModule:
  #testblock "findHomopolymerRuns":
  #  let x = findHomopolymerRuns("AACCCTTTTA")
  #  doAssert x == @[2, 1, 3, 1, 1, 4, 1, 1, 1, 1]


  echo "OK: all tests passed"