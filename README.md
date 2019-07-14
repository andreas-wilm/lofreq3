# LoFreq Version 3

Run `nimble test` to run tests (see ./tests)

Run `nimble build` to build the binary (see `./lofreq`)

# Notes for Users

LoFreq comes with one binary called `lofreq`. This binary supports sub-commands (just like samtools
does). Call the binary without any arguments for usage information, like parameters etc.

Variant calling is done in two stages:
1. Pileup
1. Calling


## Pileup

The pileup step extracts all information relevant to LoFreq from your BAM file. It basically creates a quality histogram per genome position. Run `lofreq pileup` to get usage information. All read- and base-level filtering happens at this step. The output of this step is then used for the actual variant calling return (see below).

We advise against excessive filtering, because it can bias results. Keep in mind that LoFreq was designed to model
sequencing (and mapping) errors!

The pileup also performs merging of qualities, e.g. mapping and base quality. The output file stores one quality per base. The idea for quality merging is as follows: all qualities stored in a BAM file are (or should in theory be) Phred scaled error probabilties (for example a base quality of 20 means, there is a 1% chance that this is in fact another base). So to combine mapping $$p_m$$, alignment $$p_a$$ and base qualities $$p_bq$$ you can do the following:

$$p_c = p_m + (1-p_m)*p_a + (1-p_m)*(1-p_a)*p_b$$

In other words: either this read (and therefore the base in question) is wrongly aligned ($$p_m$$) or if
it's not ($$1-p_m$$) then it might be wrongly (base-)aligned ($$p_a$$) and if it's neither mismapped or
misaligned $$(1-p_m)*(1-p_a)$$ then it might still be the wrong base ($$p_b$$).


Depending on the type of data you have you might want to ignore bases with base quality 2. This was used by Illumina for a while as so called "Read Segment Quality Control Indicator". Here, a base quality of 2 is in fact not a quality, but just means that the sequencing machine wasn't sure about the entire segment. Since it's not a base quality (but just an indicator) this affects LoFreq.


## Call

This steps takes the pileup json file as input, calls variants and outputs a VCF file. Note, no
base or read filtering will happen at this stage. Basic variant filtering is applied. See `lofreq call` for more information





# Notes for Developers

Read level filtering can only happen at the pileup stage, so it makes
sense to apply base level filtering there as well. In other words, the
variant calling routines only deal with filtered data. The data
exchanged between `plp` and `call` is basically an annotated quality
histogram per event (base or indel). The histogram allows a dense
representation for ultra high coverage data, while keeping all data
needed for variant calling. The downside is, we cannot keep multiple
qualities per vent, because the linkage is broken in a histogram. So
quality joining/mergign has to happen in the pileup phase as well. As
little filtering as possible should happen in the pileup phase,
because firstly LoFreq is meant to deal with noise and secondly you
otherwise skew results (coverage etc.)


# To Do List

- pileup: allow to specify region
- pileup: allow pileup without reference
- pileup: don't parse reference for chroms without data (currently parses all in header). best to
  parse fai to algorithm.pileup
- pileup: fix `../lofreq3.git/lofreq pileup -b chr20-40000001-45000000.bam -f GRCh38_full_analysis_set_plus_decoy_hla.fa >/dev/null`
  - beginning <= position` The file is not sorted: 40002361 40002362 [AssertionError]
  - `47ea7271-d86d-4e27-b421-39b95e056f89_Basecall_Alignment_template        16      chr20    40002363        0       3078S3I...`
  - classify as invalid CIGAR?
- add tests for invalid cigar ops: ^[ID] and ^[S][ID}
