#!/bin/bash

set -euo pipefail

LOFREQ="../lofreq"
REFFA="denv2-simulation/denv2-refseq.fa"
BAM="denv2-simulation/denv2-10haplo.bam"
REG="gi|158976983|ref|NC_001474.2|:1-10723"
VCFEVAL=../vcfeval
TRUTHVCF=denv2-simulation/denv2-10haplo_true-snp.vcf.gz

odir=$(mktemp -d)
echo "Running tests in $odir"
outvcf=$odir/denv2-10haplo.vcf.gz
log=$odir/denv2-10haplo.log
/usr/bin/time -v $LOFREQ call -a 0.0 -f $REFFA -b $BAM -r "$REG" 2>$log | bgzip > $outvcf
tabix $outvcf

testvcf=$odir/denv2-10haplo.flt.vcf.gz
bcftools view -i 'QUAL>60 && SB<60' $outvcf -O z -o $testvcf
tabix $testvcf

# Lost var wrt 2.1.4:
#  gi|158976983|ref|NC_001474.2|   8755    .       A       C       20      PASS	DP=8852;AF=0.001017;SB=0;DP4=4458,4383,5,4
# Gets qual 18.

# 2.1.4 tests allowed 15/19 missing without/with BAQ
$VCFEVAL -r $TRUTHVCF -t $testvcf --maxFN 15 --maxFP 0 -v snp

echo "PASS"

