# auto-generate using c2nim  bam_md_ext.h, "_t" removed, compile pragmas and importc added 
{.compile: "bam_md_ext.c".}
{.compile: "kprobaln_ext.c".}
{.passL: "-lm"}

type
  bam_lf_t* {.bycopy.} = object
    cigar*: ptr uint32
    qual*: ptr uint8
    seq*: ptr uint8
    pos*: int32
    l_qseq*: int32
    n_cigar* {.bitsize: 16.}: uint32

type aln_qual_strgs* = object
  ai_str*: string
  ad_str*: string
  baq_str*: string

                   
proc bam_prob_realn_core_ext(b: ptr bam_lf_t; `ref`: cstring; baq_flag: cint;
                             baq_extended: cint; idaq_flag: cint; baq_str: cstring;
                             ad_str: cstring; ai_str: cstring): cint {.cdecl, importc: "bam_prob_realn_core_ext".}


proc bam_prob_realn_core_ext*(b: ptr bam_lf_t; `ref`: cstring; baq_flag: cint;
                             baq_extended: cint; idaq_flag: cint; aqs: aln_qual_strgs): int =
  result = bam_prob_realn_core_ext(b, `ref`, baq_flag, baq_extended, idaq_flag, aqs.baq_str, aqs.ad_str, aqs.ai_str)



when isMainModule:
  echo "OK: no tests"