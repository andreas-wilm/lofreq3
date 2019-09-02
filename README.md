# LoFreq Version 3


This version is a full reimplementation of [LoFreq version 2](http://csb5.github.io/lofreq/) in Nim. tHIS  version is currently still under development, but already provides some basic and useful functionlity. For example, it comes with a multipurpose pileup engine, that generates a quality histogram for mis/matches and indels per position in JSON format.

Why Nim? Because Nim has an intuitive and clean syntax, it looks similar to Python and compiles via C to small and fast binaries. And it's simply fun.
Thanks to Brent Petersen we have a [htslib library for Nim](https://github.com/brentp/hts-nim).


## Why a new version

The old code was hard to maintain and hard to change. Some rules were confusing. We also wanted a complete separation of the pileup and the variant calling process. 


## LoFreq Basics

**FIXME**

- Qualities are error probabilities. Ideally calibrated. Hard for MQ (see [Lee et al. 2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090581)
- Assumptions: positions are independent
- Preprocessing crucial
- Originally built for viral/bacterial data for which GATK best practices is not recommended

## Compilation

To compile from source, you will need to install [Nim](https://nim-lang.org/install.html) first.

Run `nimble test` to run tests (see ./tests)

Run `nimble build` to build the binary (see `./lofreq`)

To execute the compiled library, you will need the [htslib library](https://github.com/samtools/htslib) installed
and in your LD_LIBRARY_PATH as well. Easiest way is to use the bioconda channel
and:
   
    conda install htslib
    export LD_LIBRARY_PATH=YOUR_MINICONDA_PATH/lib:$LD_LIBRARY_PATH


### Dependencies

Installation of dependencies is taking care of by Nimble.

This version is know to work with:
cligen@0.9.37, hts@0.2.19, tempfile@0.1.7

# Usage

LoFreq comes with one binary called `lofreq`. This binary supports sub-commands (just like samtools
does). Call the binary without any arguments for usage information, like parameters etc.

Variant calling is done in two stages:
1. Pileup
1. Calling


## Pileup

The pileup step extracts all information relevant to LoFreq from your BAM file.
It basically creates a quality histogram per genome position. Run `lofreq
pileup` to get usage information. All read- and base-level filtering happens at
this step. The output of this step is then used for the actual variant calling
return (see below).

We advise against excessive filtering, because it can bias results. Keep in
mind that LoFreq was designed to model sequencing (and mapping) errors!

The pileup also performs merging of qualities, e.g. mapping and base quality, an idea
intrinsic to LoFreq version 2. The output file stores one quality per base. The idea for quality merging is as
follows: all qualities stored in a BAM file are (or should in theory be) Phred
scaled error probabilties (for example a base quality of 20 means, there is a
1% chance that this is in fact another base). So to combine mapping `p_m`,
alignment `p_a` and base qualities `p_bq` you can do the following:

`p_c = p_m + (1-p_m)*p_a + (1-p_m)*(1-p_a)*p_b`

In other words: either this read (and therefore the base in question) is
wrongly mapped (`p_m`) or if it's not (`1-p_m`) then the base/indel might be wrongly
aligned (`p_a`) and if it's neither mismapped or misaligned
`(1-p_m)*(1-p_a)` then it might still be the wrongly read base/indel (`p_b`).


### Base quality 2: Illumina's Read Segment Quality Control Indicator

A base quality of value 2 is used by Illumina as so called "Read Segment
Quality Control Indicator". Here, the quality of 2 is in fact not a quality or error probability,
but just means that the sequencing machine wasn't sure about an entire
segment in the read. This obviously affects LoFreq, which treat these as bases/indels with high error probabilty.


### Some more background on the design

Read level filtering can only happen at the pileup stage, so it makes
sense to apply base level filtering there as well. In other words, the
variant calling routines only deal with filtered data. The data
exchanged between `plp` and `call` is basically an annotated quality
histogram per event (base or indel). The histogram allows a dense
representation for ultra high coverage data, while keeping all data
needed for variant calling. The downside is, we cannot keep multiple
qualities per event, because the linkage is broken in a histogram. So
quality joining/merging has to happen in the pileup phase as well. As
little filtering as possible should happen in the pileup phase,
because firstly LoFreq is meant to deal with noise and secondly you
otherwise skew results (coverage etc.)

## Call

This steps takes the pileup json file as input, calls variants and outputs a
VCF file. Note, no base or read filtering will happen at this stage. Default
variant filtering and allele frequency filtering applied. See `lofreq call` for
default values. Strand bias (SB) is reported by not filtered on. Unless your
samples was highly PCR amplified, we suggest to filter on strand bias

- FIXME: bcftools command
- FIXME: paper explaining strand bias it

## Preprocessing

LoFreq depends on quality values and a raw BAM file will usually only contain base qualities. Other qualities (indel qualities, alignment qualities like BAQ) need to be added to the BAM file if you want to use them. Ideally they are calibrated as well. 

This version of LoFreq doesn't implement any preprocessing steps. Please refer to LoFreq version 2.

A recommended preprocessing workflow looks as follows:
- Align reads with BWA 
- Quality aware realignment with LoFreq2's `viterbi`
- Insertion of indel qualities with LoFreq2's `indelqual` (also done by GATK's BQSR)
- Base quality recalibration with GATK's BQSR 
- Insertion of alignment qualities with LoFreq2's `alnqual`


# Citation

If you use LoFreq, please cite the original publication:

[Wilm _et al._ LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Res. 2012; 40(22):11189-201.](https://www.ncbi.nlm.nih.gov/pubmed/23066108)

# To Do List

- CI on Github
- Tests, tests, tests
  - If nimble test is too limiting: https://github.com/ryanlayer/ssshtest
  - add tests for invalid cigar ops: `^[ID]` and `^[S][ID]`
  - test merged qualities against old lofreq
  - full test against spike in data
  - coverage and SB output in vcf
- Performance
  - Pileup slow on nanopore data: even without printing json, dequeue initial size 100000 and release mode
- Multithreaded pileup and call as one
- Documentation:
  - Note on invalid CIGAR ops
  - Explain how it works
  - Why we force region and how to generate it
- Enable merging of alignment qualities for bases, ins and dels in processor.nim 
- Enable use of indel quals in processor.nim 
- Publish on Github and create issues
- Add Snakemake workflow and container for preprocessing with old LoFreq. Consider adding BAQ and indel multiplexed before sorting
