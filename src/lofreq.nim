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
                   "minAF": "minimum variant frequency for variants (applied at calling stage)",
                   "minBQ": "ignore bases with base quality below this value (applied at pileup stage)",
                   "maxCov": "ignore positions with coverage above this value (applied at pileup stage)",
                   "minCov": "ignore positions with coverage below this value (applied at pileup stage)",
                   "minVarQual": "minimum variant quality (applied at calling stage)",
                   "noMQ": "ignore mapping quality (applied at pileup stage)",
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
