#!/bin/bash

set -euo pipefail

# This uses the non optimized version and is hence slow
LOFREQ="../lofreq"
LOFREQ="lofreq"
FQ1=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae_comb_1.fq.gz
FQ2=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae_comb_2.fq.gz
REFFA=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae.fasta


VCFEVAL=../vcfeval
TRUTHVCF=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae_comb.laln.vcf.gz

outdir=$(mktemp -d)
echo "Running tests in $outdir"

outbam=$outdir/lf3.bam
{ bwa mem -t 8 $REFFA $FQ1 $FQ2 | \
    $LOFREQ viterbi -f $REFFA -b - | \
    samtools fixmate - - | \
    samtools sort - | \
    $LOFREQ indelqual -f $REFFA -b - | \
    $LOFREQ alnqual -f $REFFA  -b - | \
    samtools view -b - -o $outbam; } >& ${outbam}.log;
samtools index $outbam

rawvcf=$outdir/lf3.vcf.gz
$LOFREQ call -f $REFFA -b $outbam | bgzip > $rawvcf
tabix $rawvcf

fltvcf=$outdir/lf3.flt.vcf.gz
# FIXME filter should be dynamic and different per type
bcftools filter -i 'QUAL>=72 && DP>=10 && SB<30 && AF>=0.005' $rawvcf -O z -o $fltvcf
tabix $fltvcf

# FIXME update as needed
$VCFEVAL -r $TRUTHVCF -t $fltvcf -v snp --minRecall 0.45 --minPrecision 0.99
$VCFEVAL -r $TRUTHVCF -t $fltvcf -v indel --minRecall 0.61 --minPrecision 0.99

echo "PASS"
rm -rf $odir
