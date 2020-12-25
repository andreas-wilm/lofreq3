## LoFreq: variant calling routines
##
## - Author: Andreas Wilm <andreas.wilm@gmail.com>
## - License: The MIT License

 
# standard
import tables
import strutils
#import strformat

# third party
import hts
#import hts/bam/cigar

# project specific
#import utils
import bam_md_ext


const AI_TAG* = "ai"
const AD_TAG* = "ad"
const BAQ_TAG* = "lb"


# from htsnim. private there
template bam_get_seq(b: untyped): untyped =
#  cast[CPtr[uint8]](cast[uint]((b).data) + uint(((b).core.n_cigar shl 2) + (b).core.l_qname))
  cast[ptr uint8](cast[uint]((b).data) + uint(((b).core.n_cigar shl 2) + (b).core.l_qname))
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
  
 
proc createRec(rec: Record, aqs: aln_qual_strgs): string =
    # We can't set BAM values in htsnim (argh), so convert to string / SAM
    # as in https://github.com/brentp/bamject/blob/master/src/bamject.nim

    when not defined(release):
      var query: string
      discard rec.sequence(query)
      if len(aqs.ai_str)>0:
        assert len(aqs.ai_str) == len(query)
      if len(aqs.ad_str)>0:
        assert len(aqs.ad_str) == len(query)

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
    if len(aqs.ai_str)>0:
      recSplit.add(AI_TAG & ":Z:" & aqs.ai_str)
    if len(aqs.ad_str)>0:
      recSplit.add(AD_TAG & ":Z:" & aqs.ad_str)
    if len(aqs.baq_str)>0:
      recSplit.add(BAQ_TAG & ":Z:" & aqs.baq_str)

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

    # No need to delete existing tags here (done in createRecords).
    # This way we could reuse existing tags in the c function
    # if needed
    
    # Load reference if not cached
    var chrom = rec.chrom
    if not refs.hasKey(chrom):
      #stderr.writeLine("DEBUG Loading " & chrom)
      refs[chrom] = fai.get(chrom)

    const baq_flag = 1
    const baq_extended = 1
    const idaq_flag = 1
    var query: string
    discard rec.sequence(query)
    

    # bam_lf_t is actually an abstraction of bam1_t that's used everywhere in htslib.
    # I couldn't reuse this here without having to link against htslib and include the
    # header and couldn't copy the htsnim definitions here because of weird
    # name clashes "required type for b: ptr bam1_t but expression 'rec.b' is of type:
    # ptr bam1_t"
    var bam_lf: bam_lf_t
    bam_lf.pos = cast[int32](rec.start)
    bam_lf.l_qseq = rec.b.core.l_qseq
    bam_lf.n_cigar = rec.b.core.n_cigar
    bam_lf.cigar = bam_get_cigar(rec.b)
    bam_lf.qual = bam_get_qual(rec.b)
    bam_lf.seq = bam_get_seq(rec.b)
    var aqs: aln_qual_strgs
    aqs.ai_str = newString(len(query))
    aqs.ad_str = newString(len(query))
    aqs.baq_str = newString(len(query))
    
    var rc = bam_prob_realn_core_ext(addr bam_lf, refs[chrom], 
                            baq_flag, baq_extended, idaq_flag, aqs)
    doAssert rc == 0

    # fix overallocated strings.
    if aqs.ad_str[0] == '\0':
      aqs.ad_str.setlen(0)
    if aqs.ai_str[0] == '\0':
      aqs.ai_str.setlen(0)
    if aqs.baq_str[0] == '\0':
      aqs.baq_str.setlen(0)
    echo createRec(rec, aqs)
    #obam.write(createRec(rec, ad_str, ai_str, baq_str))

  #oBam.close()
  
  
when isMainModule:
  #testblock "findHomopolymerRuns":
  #  let x = findHomopolymerRuns("AACCCTTTTA")
  #  doAssert x == @[2, 1, 3, 1, 1, 4, 1, 1, 1, 1]
  echo "OK: no tests"
