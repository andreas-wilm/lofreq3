import cligen

import lofreqpkg/call as lofreq_call
# cannot be called 'call', otherwise you get an error from cligen:
# "Error: formalParams requires a proc argument."
    
when isMainModule:
  dispatch_multi([call])
