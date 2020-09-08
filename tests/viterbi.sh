#!/bin/bash

#set -eux
set -eu

# a BAM with reads containing indels
inbam=NC_000912_Mpneumoniae/viterbi/NC_000912_Mpneumoniae_comb.srt.indelonly.bam
# fasta
fasta=NC_000912_Mpneumoniae/NC_000912_Mpneumoniae.fasta
# the reference BAM produced with 2.1.5 viterbi
refbam=NC_000912_Mpneumoniae/viterbi/NC_000912_Mpneumoniae_comb.srt.indelonly.origviterbi.bam
test -e $refbam
# the test bam
testbam=${refbam%origviterbi.bam}.newviterbi.bam

# run new viterbi
../lofreq viterbi -f $fasta -i $inbam | samtools view -b -o $testbam
# count differences compared to old viterbi, but ignore tags
ndiff=$(diff <(samtools view $refbam | cut -f -11) \
    <(samtools view $testbam | cut -f -11) | grep -c '^>')
nreads=$(samtools view -c $inbam)
threshold=0.002
if [ 1 -eq "$(echo "${ndiff}/${nreads} > $threshold" | bc -l)" ]; then
    echo "FAIL: old and new viterbi implementations differ more than expected"
    exit 1
else
    echo "OK: minimal differences between old and new viterbi implementations"
fi

rm $testbam