/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009 by Heng Li <lh3@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef BAM_MD_EXT_H
#define BAM_MD_EXT_H


/* when reusing bam1_t:
        ... /home/wilma/genomics/lofreq3.git/src/lofreqpkg/alnqual.nim(159, 37) Error: type mismatch: got <ptr bam1_t, string, int literal(1), int literal(1), int literal(1), string, string, string>
        ... but expected one of:
        ... proc bam_prob_realn_core_ext(b: ptr bam1_t; refsq: cstring; baq_flag: cint;
        ...                             baq_extended: cint; idaq_flag: cint; baq_str: cstring;
        ...                             ai_str: cstring; ad_str: cstring): int
        ...   first type mismatch at position: 1
        ...   required type for b: ptr bam1_t                  <---
        ...   but expression 'rec.b' is of type: ptr bam1_t    <--- ???
        ... expression: bam_prob_realn_core_ext(rec.b, refs[chrom], 1, 1, 1, baq_str, ai_str, ad_str)
   is that because I had to redefine this here in alnqual.nim (copied from (the private) hts/private/hts_concat.nim?  
*/

typedef struct bam_lf {
   uint32_t *cigar;
   uint8_t *qual;
   uint8_t *seq;
   int32_t pos;
   int32_t l_qseq;
   uint32_t n_cigar:16;
} bam_lf_t;



int bam_prob_realn_core_ext(const bam_lf_t *b, /*const int32 *qseq, const uint8 *qual, const uint32 n_cigar, const int32 l_qseq,*/ 
                            const char *ref, 
                            int baq_flag, int baq_extended,
                            int idaq_flag, 
                            char *baq_str, char *ai_str, char *ad_str);

#endif
