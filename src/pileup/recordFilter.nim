import hts

const PAIRED_READ           = 1
const MAPPED_IN_PROPER_PAIR = 2
const READ_UNMAPPED         = 4
const MATE_UNMAPPED         = 8
const READ_REVERSE_STRAND   = 16
const MATE_REVERSE_STRAND   = 32
const FIRST_IN_PAIR         = 64
const SECOND_IN_PAIR        = 128
const NOT_PRIMARY_ALIGNMENT = 256
const FAILS_VENDOR_CHECK    = 512
const PCR_OR_DUPLICATE      = 1024
const SUPPLEMENTARY         = 2048

const DEFAULT_IGNORE_FLAGS: uint16 = 
  READ_UNMAPPED or 
  NOT_PRIMARY_ALIGNMENT or
  FAILS_VENDOR_CHECK or
  PCR_OR_DUPLICATE or
  SUPPLEMENTARY

type RecordFilter = ref object
  bam: Bam
  ignoreFlag: uint16
  chromosomeName: string

proc newRecordFilter*(bam: Bam, chromosomeName: string, 
                      ignoreFlags: varargs[uint16]): RecordFilter =
  var finalFlag: uint16

  if ignoreFlags.len == 0:
    finalFlag = DEFAULT_IGNORE_FLAGS
  else:
    for flag in ignoreFlags: 
      finalFlag = finalFlag or flag
  
  return RecordFilter(bam: bam, chromosomeName: chromosomeName, 
               ignoreFlag: finalFlag)


iterator items*(self: RecordFilter) : Record = 
  for read in self.bam.querys(self.chromosomeName):
    if (read.flag and self.ignoreFlag) == 0:
      yield read