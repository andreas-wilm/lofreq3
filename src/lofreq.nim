# standard
import cligen
#import nimprof
# third party
# project specific
import lofreqpkg/call as lofreq_call
import lofreqpkg/pileup/pileup as lofreq_pileup
# cannot be called 'call', otherwise you get an error from cligen:
# "Error: formalParams requires a proc argument."

when isMainModule:
  dispatch_multi(
    [callFromPlp, help={"plpFname": "pileup file name (LoFreq JSON format)",
                        "minVarQual": "minimum variant quality",
                        "minAF": "minimum variant frequency"}],
    [pileup, help={"bamFname": "BAM file",
                   "faFname": "fasta reference (indexed)",
                   "json": "print pileup as json (don't call variants; see also 'pretty')",
                   "minAF": "minimum variant frequency for variants",
                   "minBQ": "ignore bases with base quality below this value",
                   "maxCov": "ignore positions with coverage above this value",
                   "minCov": "ignore positions with coverage below this value",
                   "minVarQual": "minimum variant quality",
                   "noMQ": "ignore mapping quality",
                   "pretty": "pretty JSON output (cannot be used with callNow)",
                   }]
  )
