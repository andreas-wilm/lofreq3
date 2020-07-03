#!/bin/bash

set -eux

# a BAM with reads containing indels
inbam=NC_000912_Mpneumoniae/viterbi/NC_000912_Mpneumoniae_comb.srt.indelonly.bam
# fasta
fasta=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae.fasta

# the reference BAM produced with 2.1.5 viterbi
refbam=NC_000912_Mpneumoniae/indelqual/NC_000912_Mpneumoniae_comb.srt.indelonly.origindelqual-uni3040.bam
test -e $refbam
testbam=${refbam%origindelqual-uni3040.bam}newindelqual-uni3040.bam
../lofreq indelqual -f $fasta -i $inbam -u 30,40 | samtools view -b -o $testbam
set +e
diff -q <(samtools view $refbam) <(samtools view $testbam)
if [ $? -ne 0 ]; then
    echo "FAIL: old and new indelqual uniform implementations differ"
    exit 1
else
    echo "OK: no differences between old and new indelqual uniform implementations"
fi
set -e
rm $testbam

# the reference BAM produced with 2.1.5 viterbi
refbam=NC_000912_Mpneumoniae/indelqual/NC_000912_Mpneumoniae_comb.srt.indelonly.origindelqual-dindel.bam
test -e $refbam
testbam=${refbam%origindelqual-dindel.bam}indelonly.newindelqual-dindel.bam
../lofreq indelqual -f $fasta -i $inbam | samtools view -b -o $testbam
set +e
diff -q <(samtools view $refbam) <(samtools view $testbam)
if [ $? -ne 0 ]; then
    echo "FAIL: old and new indelqual dindel implementations differ"
    exit 1
else
    echo "OK: no differences between old and new indelqual dindel implementations"
fi
set -e
rm $testbam

