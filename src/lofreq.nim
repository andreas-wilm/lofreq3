import cligen

import lofreqpkg/call as lofreq_call
import lofreqpkg/pileup/pileup as lofreq_pileup
# cannot be called 'call', otherwise you get an error from cligen:
# "Error: formalParams requires a proc argument."
    
when isMainModule:
  dispatch_multi(
    [call, help={"plpFname": "Pileup file name (LoFreq JSON format)"}],
    [pileup, help={"bamFname": "BAM file",
                    "faFname": "Fasta reference (indexed)"}]
  )
