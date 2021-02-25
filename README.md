# LoFreq Version 3

LoFreq Version 3 is a reimplementation of [LoFreq version
2](http://csb5.github.io/lofreq/) in [Nim](https://nim-lang.org/).

After almost 10 years the old code had acquired a lot of technical depth, which eventually limited its extensibility and also affected usage of the program. We furthermore wanted a modular design with complete separation of the pileup and the variant calling process.

We chose [Nim](https://nim-lang.org/) for a reimplementation, because Nim has an intuitive and clean syntax. It looks similar to Python and compiles via C to small and fast binaries. And it's simply fun. Thanks to Brent Petersen the required [htslib library for
Nim](https://github.com/brentp/hts-nim) exists. We hope that this reimplementation ensures that the LoFreq development continues.

This version is fully functional, however, the pileup process is currently much slower than in the old version and will be rewritten soon.

## Citation

If you use LoFreq, please cite the original publication:

[Wilm _et al._ LoFreq: A sequence-quality aware, ultra-sensitive variant caller
for uncovering cell-population heterogeneity from high-throughput sequencing
datasets. Nucleic Acids Res. 2012;
40(22):11189-201.](https://www.ncbi.nlm.nih.gov/pubmed/23066108)

## Table of contents

- [LoFreq Version 3](#lofreq-version-3)
  - [Citation](#citation)
  - [Table of contents](#table-of-contents)
  - [For the impatient](#for-the-impatient)
  - [LoFreq in detail](#lofreq-in-detail)
    - [The idea](#the-idea)
    - [Preprocessing of your BAM file](#preprocessing-of-your-bam-file)
    - [Pileup](#pileup)
      - [Notes on Illumina's Read Segment Quality Control Indicator](#notes-on-illuminas-read-segment-quality-control-indicator)
    - [Variant Calling](#variant-calling)
    - [Postprocessing of variants](#postprocessing-of-variants)
  - [Installation](#installation)
    - [Installing a binary release](#installing-a-binary-release)
    - [Compilation from source](#compilation-from-source)
  - [To Do List](#to-do-list)

## For the impatient

LoFreq is called through a single binary called `lofreq`. It's main purpose is to call variants from a preprocessed BAM file.
It's target platform is Illumina, but it does work with other sequencing data as well.

You will achieve best variant calling results by following the following preprocessing steps for your reads:

1. Align your reads with a mapper that produces mapping qualities, e.g. [BWA](http://bio-bwa.sourceforge.net/)
1. Recommended: Realign the reads with indels with the base quality aware `lofreq viterbi` (will be
   slow for long read data)
1. Sort the BAM file with `samtools sort` (`viterbi` will change the alignment position of some reads)
1. Required for indel calling: Add indel qualities (BI and BD tags) `lofreq indelqual`. Note: the
   default (dindel qualities) is Illumina specific.
1. Recommended: Add base- and indel alignment qualities ([BAQ](https://pubmed.ncbi.nlm.nih.gov/21320865/) and a LoFreq equivalent for indels) with `lofreq alnqual`

As one step the above looks as follows:

    obam=your-output.bam;
    reffa=your-ref.fasta
    fq1=your-R1.fq.gz
    fq2=your-R2.fq.gz
    bwa mem $reffa  $fq1 $fq2 | \
      samtools fixmate - - | \
      lofreq viterbi -f $reffa -b - | \
      samtools sort - | \
      lofreq indelqual -f $reffa -b - | \
      lofreq alnqual -f $reffa  -b - | \
      samtools view -b - -o $obam;


Then, use `lofreq call` to call variants in the processed BAM file. The following will call variants in BAM file `aln.bam` against reference `ref.fa` at chromosome `chr` between positions `s` to `e`:

    lofreq call -b aln.bam -f ref.fa -r chr:s-e

Ideally you immediately convert the output to compressed vcf and index the file:

    lofreq call -b aln.bam -f ref.fa -r chr:s-e | bgzip > out.vcf.gz
    tabix out.vcf.gz

Finally filter variants on base quality and for example strand-bias as needed with `bcftools`.

Use  `lofreq help` or `lofreq cmd --help` to display usage information, parameters etc.


## LoFreq in detail

### The idea

LoFreq was originally developed as a quality-aware low-frequency variant caller. The main design goal was to find rare variants caused by haplotypes in viral or bacterial populations. Because it's quality aware, it is also largely parameter free and applicable to many sequencing technologies, assuming that the qualities are actually meaningful and translate into error probabilities and are calibrated. This assumption is not entirely true for example for mapping qualities (see [Lee et al. 2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090581). Another simplifying assumption is that all qualities (base-qualities, indel-qualities, mapping-qualities and alignment qualities) are independent. LoFreq merges all these qualities into one error probability and predict variants based on a poission-binomial distribution. This allows us to assign meaningful qualities to variants, which actually translate into the probability that a variant was called by chance.


### Preprocessing of your BAM file

Read mappers can make tiny mistakes on the base level, which can result in spurious variant calls. To avoid this `lofreq viterbi` for a base-quality aware realignment of reads with indels.

LoFreq makes heavy use of quality values and raw BAM files usually only contain raw base qualities and mapping qualities. We therefore recommend that you calibrate mapping qualities (e.g. GATK's BQSR, but be aware of assumption it makes) and add other qualities (indel qualities and alignment qualities) with `lofreq indelqual` and `lofreq alnqual`.

While the variant calling step itself is sequencing-technology agnostic, the three pre-processing steps above are optimized for Illumina reads and cannot be easily applied to e.g. Nanopore data.

### Pileup

The pileup step extracts all information relevant for the variant calling process from your BAM file. It basically creates an quality histogram per position in the genome. The command to create a pileup is `lofreq
call -p`. Run it with  `--help` to get usage information.

The pileup is a quality histogram per position in JSON format, which makes it directly usable by other programs. Please note that LoFreq applies quality merging (see below), so you will so only one quality per event.

All read-level filtering happens at this step. We advise against excessive filtering, because it can bias results. Keep in
mind that LoFreq was designed to model and deal with sequencing (and mapping) errors!

The pileup also performs merging of qualities, e.g. mapping and base quality, an
idea intrinsic to LoFreq version 2. The output file stores one quality per base.
The idea for quality merging is as follows: all qualities stored in a BAM file
are (or should in theory be) [Phred scaled](https://en.wikipedia.org/wiki/Phred_quality_score) error probabilities.
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
entire segment in the read. This affects LoFreq, which treats qualities as error probabilities, and thus bases with low quality as
bases with high error probability.

To deal with this problem all bases with BQ=2 are marked by the pileup function with a quality of -1 (argument `--minBQ`). And matches/mismatches with quality <0 are ignored by the variant calling routine. This effectively filters bases with BQ2, but still keeps them visible in the pileup (with quality -1).

Variant calls at position with lots of bases with BQ2 are to be taken with a grain of salt. You can get very different results depending on whether you include them (SNV call less likely because of perceived high chances of error) or not (SNV call more likely).

We discourage users from changing the default of `--minBQ 3`. LoFreq is build to deal with errors and excessive filtering will bias results.

### Variant Calling

This step calls variants and outputs a
VCF file. It is implemented in the `lofreq call` command. Default
variant quality filtering (`--minVarQual`) and allele frequency filters (`--minAF`) are applied. You can in addition filter bases below a minimum base quality (`--minBQ`)  and variants within a coverage range (`--minCov` and `--maxCov`).

We do not recommend to change default filters, unless you know exactly what you are doing. Especially `--minBQ` is often misused. Remember that LoFreq builds error probabilities into its calling model and excessive filtering will create unwanted biases.

See `lofreq call --help` for all supported parameters and default values.

Strand bias (SB) is reported by not used for filtering by default. Note that strand bias doesn't mean that one strand has more bases then the other, but that the distribution of alt and ref bases between forward and reverse strand is skewed. This is tested with Fisher's Exact test as also done in samtools.

Unless your samples were highly PCR amplified, we suggest to filter on strand bias (with e.g. bcftools).


### Postprocessing of variants

Our recommendation is to filter on quality, coverage, strand bias and possibly allele frequency.

LoFreq's quality is a phred-scaled error probabilty that the reported variant was predicted wrongly.
A quality of 20 therefore means, that there is 1% chance of error, 30 means 0.1%
etc. The more positions you test, the more likely you are to come
across false positives. Therefore you might want to apply multiple testing correction. For this you could take
the initial signifiance threshold of 0.05, divide it by the genome length (see also the fai index file of your reference sequence) and convert the resulting value into a
phred quality score. A good rule of thumb is to use a value of 60 for bacteria and a value of 90 for
a human-exome sized genome.

High strand-bias is usually a sign of false positive variants. Note that strand bias is not the imbalance
between base counts on forward and reverse strand, but the imbalance between ref and alt base counts on
forward and reverse strand. High strand-bias happens especially in highly PCRed data and even
more so in vicinity of primer sites.

Allele frequency filtering is by right not necessary, but a minimal value of 0.5% is a good default,
because LoFreq's resolution was experimentally tested up to that value.

A good multi-purpose filter command is the following:

     bcftools filter -i 'QUAL>=60 && DP>=10 && SB<30 && AF>=0.005' in.vcf.gz -o out.vcf.gz

This only keeps variants with

- quality higher or equal 60
- coverage higher or qqual 10
- strand bias lower than 30
- allele frequency higher or equal than 0.5%

There are other filters you might employ e.g. remove highly repetitive regions with
[mdust](https://github.com/lh3/mdust).

## Installation

### Installing a binary release

**TODO**

### Compilation from source

To compile from source, you will need to install [Nim](https://nim-lang.org/install.html) first.

Run `nimble build` to build the binary (see `./lofreq`)

Run `nimble test` to run tests (some of which depend on the successful build).

More extensive and longer running tests can be found in the shell files in the `./tests/` directory. 

To execute the compiled binary, you will need the [htslib
library](https://github.com/samtools/htslib) installed and in your
LD_LIBRARY_PATH as well.

    export LD_LIBRARY_PATH=YOUR_MINICONDA_PATH/lib:$LD_LIBRARY_PATH

Easiest way is to install htslib is to use [bioconda](https://bioconda.github.io/):

    conda install htslib

Installation of other dependencies is taking care of by Nimble.


