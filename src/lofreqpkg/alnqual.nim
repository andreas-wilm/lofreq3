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
    if skipRead(rec):
      #obam.write(rec)
      echo $rec.tostring()
      continue

    #if has_ins or has_del:
    # FIXME then what?
    
    # FIXME delete existing tags
    
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

#[ FIXME long comment
type
  bam_lf_t* {.bycopy.} = object
    cigar*: ptr uint32
    qual*: ptr uint8
    seq*: ptr uint8
    pos*: int32
    l_qseq*: int32
    n_cigar* {.bitsize: 16.}: uint32

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t bin:16, qual:8, l_qname:8;
    uint32_t flag:16, n_cigar:16;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;

typedef struct {
    bam1_core_t core;
    int l_data, m_data;
    uint8_t *data;
#ifndef BAM_NO_ID
    uint64_t id;
#endif
} bam1_t;

  cigar
  qual 
  and
  seq
  are the challenging ones because they are extracted from data

  #define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
  #define bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
  #define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))


bam_get_seq only used twice
  uint8_t *s, *r, *q, *seq = bam_get_seq(b), *bq;
  for (i = 0; i < c->l_qseq; ++i) s[i] = seq_nt16_int[bam_seqi(seq, i)];

  ins_seq[j] = seq_nt16_str[bam_seqi(bam_get_seq(b), y)];

  seq_nt16_table converts from ascii to int 0-15
  seq_nt16_str converts 0-15 int to IUPAC table
  seq_nt16_int converts 0-15 int to ACGTN int encoded
 

 bam_get_qual: qual later handed to kpa_ext_glocal so keep as-is

 bam_get_ciar
]#
  

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
    notImplementedError

  #oBam.close()
  
when isMainModule:
  #testblock "findHomopolymerRuns":
  #  let x = findHomopolymerRuns("AACCCTTTTA")
  #  doAssert x == @[2, 1, 3, 1, 1, 4, 1, 1, 1, 1]


  echo "OK: all tests passed"