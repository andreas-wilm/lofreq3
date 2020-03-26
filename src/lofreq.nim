# standard
import cligen
when compileOption("profiler"):
  import nimprof

# third party
# /

# project specific
# following imports need to be named otherwise cligen fails!?
import lofreqpkg/call as lofreq_call
import lofreqpkg/pileup/pileup as lofreq_pileup
import lofreqpkg/pileup/algorithm# because const in pileup.call() exposed here are defined there
import lofreqpkg/viterbi as lofreq_viterbi

when isMainModule:
  dispatch_multi(
    [viterbi,
      help = {"faFname": "fasta reference (indexed)",
              "bamIn": "BAM input (\"-\" for stdin)",
              "bamOut": "BAM output (\"-\" for stdout)"},
      short = {"faFname": 'f',
               "bamIn": 'i',
               "bamOut": 'o'}],
    [call_from_plp,
      help = {"plpFname": "pileup file name (LoFreq JSON format)",
              "minVarQual": "minimum variant quality",
              "minAF": "minimum variant frequency"}],
    [call,
      help = {"bamFname": "BAM file",
              "faFname": "fasta reference (indexed)",
              "regions": "Regions in the form of sq:s-e. Separate multiple regions with comma.",
              "bedFname": "BED file listing regions",
              "minVarQual": "minimum variant quality (applied at calling stage)",
              "minAF": "minimum variant frequency for variants (applied at calling stage)",
              "maxCov": "ignore positions with coverage above this value (applied at pileup stage)",
              "minCov": "ignore positions with coverage below this value (applied at pileup stage)",
              "minBQ": "ignore bases with base quality below this value (applied at pileup stage)",
              "noMQ": "ignore mapping quality (applied at pileup stage)",
              "pileup": "Don't call variants, but print pileup as JSON instead. See also 'pretty')",
              "pretty": "pretty JSON output (cannot be used with callNow)"},
      short = {"bamFname": 'b',
               "faFname": 'f',
               "regions": 'r',
               "bedFname": 'l',
               "pileup": 'p',
               "minAF": 'a',
               "minBQ": 'q',
               "maxCov": 'C',
               "minCov": 'c',
               "minVarQual": 'v',
               "noMQ": 'M',
               "pretty": 'P',}]
  )
