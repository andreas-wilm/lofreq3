## Implements a record filter for BAM files. The 'RecordFilter' object is a
## decorator for records received by querying the BAM file. It lets through
## all records which belong to a specified chromosome and are not marked with
## certain flags.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License


import hts

## All the possible flags for the records. Use this instead of raw numbers
# FIXME: an enum might be a better fit for this
const PAIRED_READ*           = 1
const MAPPED_IN_PROPER_PAIR* = 2
const READ_UNMAPPED*         = 4
const MATE_UNMAPPED*         = 8
const READ_REVERSE_STRAND*   = 16
const MATE_REVERSE_STRAND*   = 32
const FIRST_IN_PAIR*         = 64
const SECOND_IN_PAIR*        = 128
const NOT_PRIMARY_ALIGNMENT* = 256
const FAILS_VENDOR_CHECK*    = 512
const PCR_OR_DUPLICATE*      = 1024
const SUPPLEMENTARY*         = 2048

## All the flags which are ignored by default. If no flags are specified in the
## constructor, use this. 
const DEFAULT_IGNORE_FLAGS*: uint16 =
  READ_UNMAPPED or 
  NOT_PRIMARY_ALIGNMENT or
  FAILS_VENDOR_CHECK or
  PCR_OR_DUPLICATE or
  SUPPLEMENTARY

type RecordFilter = ref object
  ## The 'RecordFilter' object.
  bam: Bam
  ignoreFlag: uint16
  chromosomeName: string

proc newRecordFilter*(bam: Bam, chromosomeName: string,
                      ignoreFlags: varargs[uint16]): RecordFilter =
  ## Constructs a new 'RecordFilter' object. The reads that are allowed through
  ## the filter: 
  ## 1. Belong to the provided chromosome
  ## 2. Are marked with non of the specified flags (if no flag is specified,
  ## the default is used)
  var finalFlag: uint16

  if ignoreFlags.len == 0:
    finalFlag = DEFAULT_IGNORE_FLAGS
  else:
    for flag in ignoreFlags:
      finalFlag = finalFlag or flag
  
  return RecordFilter(bam: bam, chromosomeName: chromosomeName,
               ignoreFlag: finalFlag)


iterator items*(self: RecordFilter) : Record =
  ## Enables transparent iteration in for..in loops. Makes any 'RecordFilter'
  ## object an iterable. This method should in most cases be called implicitly.
  for read in self.bam.query(self.chromosomeName):
    if (read.flag and self.ignoreFlag) == 0:
      yield read
