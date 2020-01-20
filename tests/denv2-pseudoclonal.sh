#!/bin/bash

set -euo pipefail

LOFREQ="../lofreq"
REFFA="denv2-pseudoclonal/denv2-pseudoclonal_cons.fa"
BAM="denv2-pseudoclonal/denv2-pseudoclonal.bam"
BED="denv2-pseudoclonal/denv2-pseudoclonal_incl.bed"

VCFEVAL=../vcfeval
TRUTHVCF=denv2-pseudoclonal/denv2-pseudoclonal_true-snp.vcf.gz

odir=$(mktemp -d)
echo "Running tests in $odir"
outvcf=$odir/denv2-pseudoclonal.vcf.gz
log=$odir/denv2-pseudoclonal.log
/usr/bin/time -v $LOFREQ call -a 0.001 -f $REFFA -b $BAM -l "$BED" 2>$log | bgzip > $outvcf
tabix $outvcf

testvcf=$odir/denv2-pseudoclonal.flt.vcf.gz
bcftools view -i 'QUAL>60 && SB<60' $outvcf -O z -o $testvcf
tabix $testvcf
$VCFEVAL -r $TRUTHVCF -t $testvcf -v snp --maxFP 0 --minTP 236
# original test had --minTP 229

echo "PASS"

