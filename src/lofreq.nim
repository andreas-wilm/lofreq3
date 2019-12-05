# standard
import cligen
#import nimprof

# third party
# /

# project specific
import lofreqpkg/call as lofreq_call
# cannot be called 'call', otherwise you get an error from cligen:
# "Error: formalParams requires a proc argument."
import lofreqpkg/pileup/pileup as lofreq_pileup


when isMainModule:
  dispatch_multi(
    [call_from_plp, help = {"plpFname": "pileup file name (LoFreq JSON format)",
                          "minVarQual": "minimum variant quality",
                          "minAF": "minimum variant frequency"}],
    [call, help = {"bamFname": "BAM file",
                   "faFname": "fasta reference (indexed)",
                   "pileup": "Don't call variants, but print pileup as JSON instead. See also 'pretty')",
                   "minAF": "minimum variant frequency for variants",
                   "minBQ": "ignore bases with base quality below this value",
                   "maxCov": "ignore positions with coverage above this value",
                   "minCov": "ignore positions with coverage below this value",
                   "minVarQual": "minimum variant quality",
                   "noMQ": "ignore mapping quality",
                   "pretty": "pretty JSON output (cannot be used with callNow)",
                   },
           short = {"bamFname": 'b',
                    "faFname": 'f',
                    "pileup": 'p',
                    "minAF": 'a',
                    "minBQ": 'q',
                    "maxCov": 'C',
                    "minCov": 'c',
                    "minVarQual": 'v',
                    "noMQ": 'M',
                    "pretty": 'P',}]
  )
