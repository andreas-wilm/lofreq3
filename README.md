# LoFreq Version 3

Run `nimble test` to run tests (see ./tests)

Run `nimble build` to build the binary (see `./lofreq`)


# Notes for Developers

Read level filtering can only happen at the pileup stage, so it makes
sense to apply base level filtering there as well. In other words, the
variant calling routines only deal with filtered data. The data
exchanged between `plp` and `call` is basically an annotated quality
histogram per event (base or indel). The histogram allows a dense
representation for ultra high coverage data, while keeping all data
needed for variant calling. The downside is, we cannot keep multiple
qualities per vent, because the linkage is broken in a histogram. So
quality joining/mergign has to happen in the pileup phase as well. As
little filtering as possible should happen in the pileup phase,
because firstly LoFreq is meant to deal with noise and secondly you
otherwise skew results (coverage etc.)
