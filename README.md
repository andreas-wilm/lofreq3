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
  - [Usage](#usage)
    - [Pileup](#pileup)
    - [Call](#call)
  - [LoFreq explained](#lofreq-explained)
    - [The idea](#the-idea)
    - [Overview of the workflow](#overview-of-the-workflow)
    - [Preprocessing of your BAM file](#preprocessing-of-your-bam-file)
    - [Pileup](#pileup-1)
      - [Note on Illumina's Read Segment Quality Control Indicator](#note-on-illuminas-read-segment-quality-control-indicator)
    - [Variant Calling](#variant-calling)
    - [Postprocessing of variants](#postprocessing-of-variants)
  - [Installation](#installation)
    - [Installing a binary release](#installing-a-binary-release)
    - [Using a container](#using-a-container)
    - [Compilation from source](#compilation-from-source)
- [To Do List](#to-do-list)

## For the impatient

**FIXME**

## Usage

LoFreq comes with one binary called `lofreq`. This binary supports sub-commands, currently `call` and `pileup`. Run `lofreq help` or `lofreq cmd --help` to display usage information, parameters etc.

Variant calling is done in two stages: first you generate a pileup and that's input to the calling routines.

### Pileup

**FIXME**

### Call

**FIXME**

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
6. [Pileup with `lofreq pileup`](#pileup)
7. [Variant calling with `lofreq call`](#calling)
8. [Postprocessing](#postprocessing-of-variants)

### Preprocessing of your BAM file

Steps prior to the pileup stage are usually referred to as preprocessing.

LoFreq depends on quality values and a raw BAM file will usually only contain
raw base qualities and mapping qualities. Ideally you calibrate mapping qualities and add other qualities (indel qualities, alignment qualities like BAQ) to the BAM file.

This version of LoFreq doesn't implement any preprocessing steps. Please refer
to LoFreq version 2.

**FIXME**

### Pileup

The pileup step extracts all information relevant for the variant calling process from your BAM file. It basically creates an quality histogram per position in the genome. The command to create a pileup is `lofreq
pileup`. Run it with  `--help` to get usage information.

All read-level filtering happens at
this step. We advise against excessive filtering, because it can bias results. Keep in
mind that LoFreq was designed to model and deal with sequencing (and mapping) errors!

The pileup also performs merging of qualities, e.g. mapping and base quality, an
idea intrinsic to LoFreq version 2. The output file stores one quality per base.
The idea for quality merging is as follows: all qualities stored in a BAM file
are (or should in theory be) [Phred scaled](FIXME) error probabilties. For example, a base
quality of 20 means, there is a 1% chance that this is in fact another base. To combine mapping `p_m`, alignment `p_a` and base qualities `p_bq` you can do
the following:

`p_c = p_m + (1-p_m)*p_a + (1-p_m)*(1-p_a)*p_b`

In other words: either this read (and therefore the base in question) is
wrongly mapped (`p_m`) or if that's not the case (`1-p_m`) then the base/indel might be wrongly
aligned (`p_a`) and if it's neither mismapped nor misaligned
`(1-p_m)*(1-p_a)` then it might still be a wrongly read base/indel (`p_b`).

#### Note on Illumina's Read Segment Quality Control Indicator

A base quality of value 2 is used by Illumina as so called "Read Segment Quality
Control Indicator". Here, the quality of 2 is in fact not a quality or error
probability, but just means that the sequencing machine wasn't sure about an
entire segment in the read. This obviously affects LoFreq, which treat these as
bases/indels with high error probabilty.

### Variant Calling

This steps takes the pileup json file as input, calls variants and outputs a
VCF file. It is implemented in the `lofreq call` command. Default
variant filtering and allele frequency filters are applied. See `lofreq call --help` for default values.

Strand bias (SB) is reported by not used for filtering by default. Note that strand bias doesn't mean that one strand has more bases then the other, but that the distribution of alt anf ref bases between forward and reverse strand is skewed. This is tested with Fisher's Exact test as also done in samtools.

Unless your samples were highly PCR amplified, we suggest to filter on strand bias.

- FIXME: bcftools command
- FIXME: paper explaining strand bias it

### Postprocessing of variants

**FIXME**: idea only

## Installation

### Installing a binary release

**FIXME**

### Using a container

### Compilation from source

To compile from source, you will need to install [Nim](https://nim-lang.org/install.html) first.

Run `nimble build` to build the binary (see `./lofreq`)

Run `nimble test` to run tests (some of which depend on the successful build).

To execute the compiled library, you will need the [htslib
library](https://github.com/samtools/htslib) installed and in your
LD_LIBRARY_PATH as well.

    export LD_LIBRARY_PATH=YOUR_MINICONDA_PATH/lib:$LD_LIBRARY_PATH

Easiest way is to install htslib is to use [bioconda](FIXME):

    conda install htslib

Installation of other dependencies is taking care of by Nimble.

This version is know to work with:
cligen@0.9.37, hts@0.2.19, tempfile@0.1.7

# To Do List

In order of importance:

- testing: base qualities, ins qualities, del qualities
- CI on Github
- Wrapper: Multithreaded pileup and call as one, with region inference
- Create release and announce availability
- Tests
  - If nimble test is too limiting: https://github.com/ryanlayer/ssshtest
  - add tests for invalid cigar ops: `^[ID]` and `^[S][ID]`
  - test merged qualities against old lofreq
  - full test against spike in data
  - coverage and SB output in vcf
- Add Snakemake workflow
- Add workflow and container for end to end processing with old LoFreq. Consider adding BAQ and indel multiplexed before sorting
- Performance
  - Pileup slow on nanopore data: even without printing json, dequeue initial size 100000 and release mode
- Extend documentation:
  - Note on invalid CIGAR ops
  - Explain how it works in detail
  - Why we force region and how to generate it
  - FAQ
  - Split usage from explanation and add section for the impatient
