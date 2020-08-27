# Code

- Branch: alnqual
- Source file: ./src/lofreqpkg/alnqual.nim

# Example run

    fa=../lofreq3-testdata/NC_000912_Mpneumoniae/NC_000912_Mpneumoniae.fasta
    bam=../lofreq3-testdata/NC_000912_Mpneumoniae/NC_000912_Mpneumoniae_comb.srt.bam
    ./lofreq alnqual -f $fa -i $bam


# Previous implementation

- lofreq_alnqual.c:main_alnqual() calls bam_md_ext.c:bam_prob_realn_core_ext() per record
- bam_md_ext.c:bam_prob_realn_core_ext() computes BAQ and IDAQ:
  1. xe, xb, bw are for start, end and bandwidth
  1. kpa_ext_glocal(r, xe-xb, s, c->l_qseq, qual, &conf, state, q, pd, &bw);
  1. BAQ computed in function
  1. IDAQ computed by calling bam_md_ext.c:idaq(bam1_t *b, const char *ref, double **pd, int xe, int xb, int bw);
- We had planned to outsource BAQ to `samtools calmd`! In that case we are computing pd again in kpa_ext_glocal() for alnqual
- Might be easiest to wrap bam_md_ext.c:bam_prob_realn_core_ext()
  - Return tags as strings (no BAM modification)
  - Should work as long as we can get access to bam1_t *b:
        # https://github.com/brentp/hts-nim/blob/master/src/hts/bam.nim
        type Record* = ref object of RootObj
        ## Record is a single alignment object.
        b*: ptr bam1_t
        hdr: Header


# Progress

## 21082020
Copied orig C files and added pragma to alnqula.nim
Commented
//#include "samutils.h"
//#include "defaults.h"
to remove further deps, but need to import needed functions and constants from there

## 22082020
Defined bam1_t in alnqual.nim by copying defintion from htsnim
ugly and results in following error anyway:
        ...   required type for b: ptr bam1_t
        ...   but expression 'rec.b' is of type: ptr bam1_t
Hand down qual, seq and cigar instead of bam1_t? needs modification all over the place, but
that awy we might be able to remove hts dep
Hand down values and then build fake struct in c to hand down to idaq(). problem are all these hts function e.g. to get seq
Leave alnqual as lofreq2 module?

## 30082020

- Imported all htslib.h/sam.h constants and defines into bam_md_ext.h
- Created bam_md_ext.nim with bam_lf as replacement of bam1_core
- bam1_core_t and bam1_t still contained in bam_md_ext.c to avoid compiler errors for now (but actual vars are NULL)


TODO
- use bam_lf whereever bam1_core_t and bam1_t where needed
- Replace functions bam_get_seq, bam_get_qual, bam_get_cigar by using bam_lf
- Remove above functions and orphaned bam1_core_t and bam1_t types from bam_md_ext.c
- Return all three tags and don't modify
