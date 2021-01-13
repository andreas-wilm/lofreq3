#!/bin/bash

set -u
#set -ux

# should really write a nim program to test this
# rather than running lofreq three times
# in bash

# a BAM with reads containing indels
inbam=NC_000912_Mpneumoniae/viterbi/NC_000912_Mpneumoniae_comb.srt.indelonly.bam
test -e $inbam
# fasta
fasta=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae.fasta
test -e $fasta
# the reference BAM produced with 2.1.5 alnqual
refbam=NC_000912_Mpneumoniae/alnqual/NC_000912_Mpneumoniae_comb.srt.indelonly.alnqual.bam
test -e $refbam
# the test bam
testbam=${refbam%alnqual.bam}.newAlnqual.bam

# run new alnqual
../lofreq alnqual -f $fasta -b $inbam | samtools view -b -o $testbam

ndiff=$(diff <(../lofreq alnqual -f $fasta -b $inbam | grep -Pow 'ai:Z:.*[^\t]' | cut -f 1) <(samtools view $refbam | grep -Pow 'ai:Z:.*[^\t]' | cut -f 1) | grep -c '^>')
if [ $ndiff -gt 0 ]; then
    echo "FAIL: old and new alnqual AI implementations differ"
    exit 1
fi

ndiff=$(diff <(../lofreq alnqual -f $fasta -b $inbam | grep -Pow 'ad:Z:.*[^\t]' | cut -f 1) <(samtools view $refbam | grep -Pow 'ad:Z:.*[^\t]' | cut -f 1) | grep -c '^>')
if [ $ndiff -gt 0 ]; then
    echo "FAIL: old and new alnqual AD implementations differ"
    exit 1
fi

ndiff=$(diff <(../lofreq alnqual -f $fasta -b $inbam | grep -Pow 'lb:Z:.*[^\t]' | cut -f 1) <(samtools view $refbam | grep -Pow 'lb:Z:.*[^\t]' | cut -f 1) | grep -c '^>')
if [ $ndiff -gt 0 ]; then
    echo "FAIL: old and new alnqual BAQ implementations differ"
    exit 1
fi

echo "OK: alnqual implementations give identical results"

rm $testbam
