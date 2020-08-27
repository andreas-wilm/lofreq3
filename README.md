# LoFreq Version 3

This version is a reimplementation of [LoFreq version
2](http://csb5.github.io/lofreq/) in [Nim](https://nim-lang.org/). The current version is still under development, but already provides basic functionlity. For
example, it comes with a multipurpose pileup engine, that generates a quality
histogram for mis/matches and indels per position in JSON format and calls variants from this input. It is however missing preprocessing modules (e.g. realignment and insertion of alignment quality) that version 2 provided.

## Why a new version

The old code had aquired a lot of technical depth, which eventually limited its extensibility and also affected usage of the program. We furthermore wanted a modular design with complete separation of the pileup and the variant calling process.

We chose [Nim](https://nim-lang.org/) for a reimplementation, because Nim has an intuitive and clean syntax. It looks similar to Python and compiles via C to small and fast binaries. And it's simply fun. Thanks to Brent Petersen the required [htslib library for
Nim](https://github.com/brentp/hts-nim) exists. We hope that this reimplementation ensures that the LoFreq development continues.

## Citation

If you use LoFreq, please cite the original publication:

[Wilm _et al._ LoFreq: A sequence-quality aware, ultra-sensitive variant caller
for uncovering cell-population heterogeneity from high-throughput sequencing
datasets. Nucleic Acids Res. 2012;
40(22):11189-201.](https://www.ncbi.nlm.nih.gov/pubmed/23066108)

## Table of contents

- [LoFreq Version 3](#lofreq-version-3)
  - [Why a new version](#why-a-new-version)
  - [Citation](#citation)
  - [Table of contents](#table-of-contents)
  - [For the impatient](#for-the-impatient)
  - [LoFreq explained](#lofreq-explained)
    - [The idea](#the-idea)
    - [Overview of the workflow](#overview-of-the-workflow)
    - [Preprocessing of your BAM file](#preprocessing-of-your-bam-file)
    - [Pileup](#pileup)
      - [Notes on Illumina's Read Segment Quality Control Indicator](#notes-on-illuminas-read-segment-quality-control-indicator)
    - [Variant Calling](#variant-calling)
    - [Postprocessing of variants](#postprocessing-of-variants)
  - [Installation](#installation)
    - [Installing a binary release](#installing-a-binary-release)
    - [Using a container](#using-a-container)
    - [Compilation from source](#compilation-from-source)
- [To Do List](#to-do-list)
  - [Testing](#testing)
  - [Documentation:](#documentation)
  - [Performance](#performance)
  - [Features](#features)
  - [Release](#release)

## For the impatient


LoFreq comes with one binary called `lofreq`. This binary supports sub-commands, currently `call` and `call_from_plp`. Run `lofreq help` or `lofreq cmd --help` to display usage information, parameters etc.

Use the following to call variants in BAM file `aln.bam` against reference `ref.fa` at chromosome `chr` between positions `s` to `e`:

    lofreq call -b aln.bam -f ref.fa -r chr:s-e


To create a pileup in JSON format (with merged qualities) run as above, but with added `-p`

    lofreq call -b aln.bam -f ref.fa -r chr:s-e -p



## LoFreq explained

**FIXME**

### The idea

**FIXME**

- Originally design to be a technology indepedent, quality-aware low-frequency variant caller. Can be used to find rare variants caused by haplotypes in viral or bacterial population.
- Because quality aware, largely parameter free and applicable to many sequencing technologies
- All sequencing related qualities are error probabilities.
- Ideally calibrated.
- Hard for MQ (see [Lee et al. 2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090581)
- One other assumption is that all positions are independent.
- Preprocessing crucial
- Originally built for viral/bacterial data for which GATK best practices is not recommended

### Overview of the workflow

A recommended preprocessing workflow looks as follows:

1. Align reads with BWA
2. Optional: Quality aware base realignment with LoFreq2's `viterbi`
3. Required for indel calling: Insertion of indel qualities with LoFreq2's `indelqual` (also done by GATK's BQSR).
4. Recommended: Base quality recalibration with GATK's BQSR
5. Recommended: Insertion of alignment qualities with LoFreq2's `alnqual`
6. [Variant calling with `lofreq call`](#calling)
7. [Postprocessing](#postprocessing-of-variants)

### Preprocessing of your BAM file

Steps prior to the pileup stage are usually referred to as preprocessing.

LoFreq depends on quality values and a raw BAM file will usually only contain
raw base qualities and mapping qualities. Ideally you calibrate mapping qualities and add other qualities (indel qualities, alignment qualities like BAQ) to the BAM file.

This version of LoFreq doesn't implement any preprocessing steps. Please refer
to LoFreq version 2.

**FIXME**

### Pileup

The pileup step extracts all information relevant for the variant calling process from your BAM file. It basically creates an quality histogram per position in the genome. The command to create a pileup is `lofreq
call -p`. Run it with  `--help` to get usage information.

The pileup is a quality histogram per position in JSON format, which makes it directly usable by other programs. Please note that LoFreq applies quality merging (see below), so you will so only one quality per event.

All read-level filtering happens at this step. We advise against excessive filtering, because it can bias results. Keep in
mind that LoFreq was designed to model and deal with sequencing (and mapping) errors!

The pileup also performs merging of qualities, e.g. mapping and base quality, an
idea intrinsic to LoFreq version 2. The output file stores one quality per base.
The idea for quality merging is as follows: all qualities stored in a BAM file
are (or should in theory be) [Phred scaled](https://en.wikipedia.org/wiki/Phred_quality_score) error probabilties.
For example, a base quality of 20 means, there is a 1% chance that this is in fact another base.
To combine mapping `p_m`, alignment `p_a` and base qualities `p_bq` you can do the following:

`p_c = p_m + (1-p_m)*p_a + (1-p_m)*(1-p_a)*p_b`

In other words: either this read (and therefore the base in question) is
wrongly mapped (`p_m`) or if that's not the case (`1-p_m`) then the base/indel might be wrongly
aligned (`p_a`) and if it's neither mismapped nor misaligned
`(1-p_m)*(1-p_a)` then it might still be a wrongly read base/indel (`p_b`).

#### Notes on Illumina's Read Segment Quality Control Indicator

A base quality (BQ) of value 2 (ASCII: `#`) is used by Illumina as so called "Read Segment Quality
Control Indicator". Here, the quality of 2 is in fact not a quality or error
probability, but just means that the sequencing machine wasn't sure about an
entire segment in the read. This affects LoFreq, which treats qualities as error probabilties, and thus bases with low quality as
bases with high error probabilty.

To deal with this problem all bases with BQ=2 are marked by the pileup function with a quality of -1 (argument `--minBQ`). And matches/mismatches with quality <0 are ignored by the variant calling routine. This effectively filters bases with BQ2, but still keeps them visible in the pileup (with quality -1).

Variant calls at position with lots of bases with BQ2 are to be taken with a grain of salt. You can get very different results depending on whether you include them (SNV call less likely because of perceived high chances of error) or not (SNV call more likely).

We discourage users from changing the default of `--minBQ 3`. LoFreq is build to deal with erros and excessive filtering will bias results.

### Variant Calling

This step calls variants and outputs a
VCF file. It is implemented in the `lofreq call` command. Default
variant quality filtering (`--minVarQual`) and allele frequency filters (`--minAF`) are applied.
See `lofreq call --help` for default values.

Strand bias (SB) is reported by not used for filtering by default. Note that strand bias doesn't mean that one strand has more bases then the other, but that the distribution of alt anf ref bases between forward and reverse strand is skewed. This is tested with Fisher's Exact test as also done in samtools.

Unless your samples were highly PCR amplified, we suggest to filter on strand bias.


### Postprocessing of variants


- **FIXME**: bcftools command
- **FIXME**: paper explaining strand bias it


## Installation

### Installing a binary release

**FIXME**

### Using a container

### Compilation from source

To compile from source, you will need to install [Nim](https://nim-lang.org/install.html) first.

Run `nimble build` to build the binary (see `./lofreq`)

Run `nimble test` to run tests (some of which depend on the successful build).

More extensive and longer running tests can be found in the shell files in the `./tests/` directory. 

To execute the compiled binary, you will need the [htslib
library](https://github.com/samtools/htslib) installed and in your
LD_LIBRARY_PATH as well.

    export LD_LIBRARY_PATH=YOUR_MINICONDA_PATH/lib:$LD_LIBRARY_PATH

Easiest way is to install htslib is to use [bioconda](FIXME):

    conda install htslib

Installation of other dependencies is taking care of by Nimble.

This version is know to work with:
cligen@0.9.37, hts@0.2.19, tempfile@0.1.7

# To Do List


## Testing

- Write snakefile for end to end testing of each step and comparison to old LoFreq (run both)
- Test merged qualities against old lofreq
- Full test against spike in data
- mincov and maxcov filter (coverage() function) in the presence of indels
- Coverage in vcf
- SB output in vcf
- Add CI on Github
- Implement filter or doc bcftools recipes (and add docs)
- Note: If nimble test is too limiting: https://github.com/ryanlayer/ssshtest

## Documentation:

- Add FAQ

## Performance

- Rewrite pileup
- Parallelism: reader queue with async calling (sidestepping json conversion)
- even the release compiled version is a few times slower then the old version (write version that
  skips that parsing just for benchmark purposes?)
- Pileup slow on nanopore data: even without printing json, dequeue initial size 100000 and release mode


## Features

- Reimplement indel qual. Challenge is to remove htslib dep and compute BAQ at the same time

## Release
- Create release with static binary and announce availability
