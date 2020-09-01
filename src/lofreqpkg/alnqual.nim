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
#import hts/bam/cigar

# project specific
import utils
import bam_md_ext

const AI_TAG = "ai"
const AD_TAG = "ad"
const BAQ_TAG = "lb"


# from htsnim. private there
template bam_get_seq(b: untyped): untyped =
#  cast[CPtr[uint8]](cast[uint]((b).data) + uint(((b).core.n_cigar shl 2) + (b).core.l_qname))
  cast[ptr uint8](cast[uint]((b).data) + uint(((b).core.n_cigar shl 2) + (b).core.l_qname))


# from htsnim. private there
template bam_get_seq(b: untyped): untyped =
#  cast[CPtr[uint8]](cast[uint]((b).data) + uint(((b).core.n_cigar shl 2) + (b).core.l_qname))
  cast[ptr uint8](cast[uint]((b).data) + uint(((b).core.n_cigar shl 2) + (b).core.l_qname))


template bam_get_qual*(b: untyped): untyped =
  #cast[CPtr[uint8]](cast[uint]((b).data) + uint(uint((b).core.n_cigar shl 2) + uint((b).core.l_qname) + uint((b.core.l_qseq + 1) shr 1)))
  cast[ptr uint8](cast[uint]((b).data) + uint(uint((b).core.n_cigar shl 2) + uint((b).core.l_qname) + uint((b.core.l_qseq + 1) shr 1)))


proc skipRead(rec: Record): bool =
  if rec.flag.secondary or rec.flag.qcfail or rec.flag.dup or rec.flag.unmapped:
    return true
  return false
  

proc createRec(rec: Record, ai: string, ad: string, baq: string): string =
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
      if fieldSplit[0] in @[AI_TAG, AD_TAG, BAQ_TAG]:
        delIndices.add(i)
    for i, j in delIndices.pairs:
      recSplit.delete(j-i)

    # and add again
    recSplit.add(AI_TAG & ":Z:" & ai)
    recSplit.add(AD_TAG & ":Z:" & ad)
    recSplit.add(BAQ_TAG & ":Z:" & baq)

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

  for rec in iBam:
    if skipRead(rec):
      #obam.write(rec)
      echo $rec.tostring()
      continue

    #if has_ins or has_del:
    # FIXME then what?
    
    # no need to delete existing tags (done in create records)
    # this way we could reuse existing tags in teh c function
    # if needed
    
    # load reference if not cached
    var chrom = rec.chrom
    if not refs.hasKey(chrom):
      #stderr.writeLine("DEBUG Loading " & chrom)
      refs[chrom] = fai.get(chrom)

    const baq_flag = 1
    const baq_extended = 1
    const idaq_flag = 1
    var query: string
    discard rec.sequence(query)
    var ai_str = newString(len(query))
    var ad_str = newString(len(query))
    var baq_str = newString(len(query))

    var bam_lf: bam_lf_t
    bam_lf.pos = cast[int32](rec.start)
    bam_lf.l_qseq = rec.b.core.l_qseq# FIXME actually len seq which we don't have here yet
    bam_lf.n_cigar = rec.b.core.n_cigar# FIXME use len(cigar(rec)) instead
    bam_lf.cigar = bam_get_cigar(rec.b)#(cast[ptr uint32](((cast[int]((rec.b).data)) + cast[int]((rec.b).core.l_qname))))
    bam_lf.qual = bam_get_qual(rec.b)
    bam_lf.seq = bam_get_seq(rec.b)# our own template might work. after all bam_get_cigar works as well
    var rc = bam_prob_realn_core_ext(addr bam_lf, refs[chrom], 
                            baq_flag, baq_extended, idaq_flag, 
                            baq_str, ai_str, ad_str)
    echo "FIXME nim baq_str=" & $baq_str
    echo "FIXME nim ai_str=" & $ai_str
    echo "FIXME nim ad_str=" & $ad_str

    echo createRec(rec, ai_str, ad_str, baq_str)
    #obam.write(createRec(rec, ai_str, ad_str, baq_str))

    notImplementedError

  #oBam.close()
  
when isMainModule:
  #testblock "findHomopolymerRuns":
  #  let x = findHomopolymerRuns("AACCCTTTTA")
  #  doAssert x == @[2, 1, 3, 1, 1, 4, 1, 1, 1, 1]


  echo "OK: all tests passed"