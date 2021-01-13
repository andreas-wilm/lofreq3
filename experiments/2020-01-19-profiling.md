# 5b17e2d

## Compile with profiling arguments and run:

    nim c  -d:nimCallDepthLimit=10000 --profiler:on --stacktrace:on src/lofreq.nim
    ./src/lofreq call -f tests/denv2-simulation/denv2-refseq.fa -b tests/denv2-simulation/denv2-10haplo.bam -r "gi|158976983|ref|NC_001474.2|:1-100"

Looking at top results:

    grep -v toolchains profile_results.txt  | grep -A3 Entry | head -n 20

    Entry: 1/282 Calls: 22/1164 = 1.89% [sum: 22; 22/1164 = 1.89%]
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/processor.nim: processMatches 1126/1164 = 96.74%
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/algorithm.nim: processEvent 1126/1164 = 96.74%
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/algorithm.nim: pileup 1135/1164 = 97.51%
    --
    Entry: 2/282 Calls: 20/1164 = 1.72% [sum: 42; 42/1164 = 3.61%]
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/storage/containers/qualityHistogram.nim: add 616/1164 = 52.92%
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/storage/containers/operationData.nim: add 616/1164 = 52.92%
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/storage/containers/positionData.nim: addDeletion 206/1164 = 17.70%
    --
    Entry: 3/282 Calls: 18/1164 = 1.55% [sum: 60; 60/1164 = 5.15%]
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/storage/containers/qualityHistogram.nim: add 616/1164 = 52.92%
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/storage/containers/operationData.nim: add 616/1164 = 52.92%
      /home/wilma/genomics/lofreq3.git/src/lofreqpkg/pileup/storage/containers/positionData.nim: addMatch 210/1164 = 18.04%

- No surprises: quality histograms additions take up most time.
- However the following routines (now shown) in insQualityAt(), baseAlnQualityAt() and delAlnQualityAt (each around 1.2%) could potentially be optimized. Htsnim:tag() also shows up in results ("memoization" of tag return?). #TODO


## Benchmark on high cov data

        Command being timed: "../lofreq.git/src/lofreq/lofreq plpsummary -B -f tests/denv2-simulation/denv2-refseq.fa tests/denv2-simulation/denv2-10haplo.bam -r gi|158976983|ref|NC_001474.2|:1-2000"
        User time (seconds): 17.50

        Command being timed: "../lofreq.git/src/lofreq/lofreq call -B -f tests/denv2-simulation/denv2-refseq.fa tests/denv2-simulation/denv2-10haplo.bam -r gi|158976983|ref|NC_001474.2|:1-2000"
        User time (seconds): 20.82

        Command being timed: "lofreq call -p -f tests/denv2-simulation/denv2-refseq.fa -b tests/denv2-simulation/denv2-10haplo.bam -r gi|158976983|ref|NC_001474.2|:1-2000"
        User time (seconds): 45.18

        Command being timed: "lofreq call -f tests/denv2-simulation/denv2-refseq.fa -b tests/denv2-simulation/denv2-10haplo.bam -r gi|158976983|ref|NC_001474.2|:1-2000"
        User time (seconds): 57.50

So pileup and call, both roughly 3x slower than 2.1.4

With compiler optimization:

`nimble install --debug --passNim:"--passC:-ffast-math"`
- pileup: little effect (44 instead of 45s)
- call: drops to 50sec from 57 (#TODO average over multiple runs)

`nimble install --debug --passNim:"--d:danger"`
- pileup: still at 44 (same as fast-math above
- call: 51 sec


