BAM contains one read. Test checks for correct position of indels and their qualities


$ samtools view tests/pileup_samples/indels/alignments.bam
1       0       ref     1       60      1M2D3M2I1M      *       0       0       AAAAGCC 5555555 BI:Z:???????    BD:Z:IIIIIII

$ samtools tview -d t tests/pileup_samples/indels/alignments.bam tests/pileup_samples/indels/ref.fa
1           11        21        31        41        51        61
AACACG**CCTTAAGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
.KK.AA  .
.**.AAGC.

