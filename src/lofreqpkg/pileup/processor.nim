## This module provides an implementation of the a processor used in the pileup
## algorithm (the module 'algorithm'). It gives the user an option to define
## provide functions for calculating base qualities of matches, insertions and
## deletions respecively, depending on specific needs.  The main idea behind
## separating the processing from the algorithm (iterating) is extensibility.
## This way, the 'Processor' can be configured to calculate qualities of all
## operation in arbitrary ways. Additionaly, if a user wants to substantially
## change the logic by which the data is processed (e.g. for a simpler form of
## pileup), they can implement a new processor and pass that to the algorithm
## instead without any need to change the algorithm itself.
##
## - Author: Filip Sodic <filip.sodic@gmail.com>
## - License: The MIT License

# standard library
#import strutils
# import logging
#import math
# third party
import hts
# project specific
import ../utils

#var logger = newConsoleLogger(fmtStr = verboseFmtStr, useStderr = true)


const DEFAULT_BLANK_QUALITY = -1

const INS_ALN_QUAL_TAG = "ai"# ins alignment quality tag. defined by LoFreq v2
const DEL_ALN_QUAL_TAG = "ad"# del alignment quality tag. Defined by LoFreq v2
const BASE_ALN_QUAL_TAG = "BQ"# base alignment quality tag (BAQ)

const INS_QUAL_TAG = "BI"# ins quality tag
const DEL_QUAL_TAG = "BD"# del quality tag


## The expected type of procedures for calculating the qualities. All relevant
## information should be obtainable through the record and the index.
## NOTE: When dealing with a deletion or an insertion, the index is actually
## the index of the base TO THE LEFT (the reference index of the last mutual
## base.  Hovewer, this should not really make a difference since, in most
## cases, all deletion and insertion quality data is found in the record's
## custom fields.
type TQualityProc = proc (r: Record, i: int, useMQ: bool): int


# for optimization purposes to avoid having to parse/decode tags multiple times
type TReadQualityBuffer* = object
    baseQuals: seq[uint8]# base qualities
    baseAlnQuals: seq[uint8]# base alignment qualities
    insQuals: seq[uint8]# insertion qualities
    insAlnQuals: seq[uint8]# insertion alignment qualities
    delQuals: seq[uint8]# deletion qualities
    delAlnQuals: seq[uint8]# deletion alignment qualties


type Processor[TStorage] = ref object
  ## The 'Processor' type. Its fields are configuration options.
  storage: TStorage
  matchQualityAt: TQualityProc# FIXME remove
  deletionQualityAt: TQualityProc# FIXME remove
  insertionQualityAt: TQualityProc# FIXME remove
  readQualityBuffer: TReadQualityBuffer# of current read
  useMQ: bool
  minBQ: int# minimum base quality. everything below will be recorded as -1.


proc mergeQuals*(q_m: int, q_a: int, q_b: int): int =
  # FIXME can we do the calculations in log space?
  # q_m = mapping quality
  # q_a = alignment quality
  # q_b = base / indel quality
  let p_m = qual2prob(q_m)# mapping error
  let p_a = qual2prob(q_a)# alignment error
  let p_b = qual2prob(q_b)# base error
  let p_c = p_m + (1-p_m)*p_a + (1-p_m)*(1-p_a)*p_b
  prob2qual(p_c)


proc insQualityAt(r: Record, i: int):  Natural =
  let iq = tag[cstring](r, INS_QUAL_TAG)
  if iq.isSome:
    return decodeQual(iq.get[i])
  else:
    return high(int)


proc insQualities(r: Record, q: var seq[uint8]): seq[uint8] =
  if len(q) != r.b.core.l_qseq:
    q.set_len(r.b.core.l_qseq)
  let iq = tag[cstring](r, INS_QUAL_TAG)
  for i in 0..<int(r.b.core.l_qseq):
    if iq.isSome:
      q[i] = decodeQual(iq.get[i])
    else:
      q[i] = high(uint8) 


proc delQualityAt(r: Record, i: int): Natural =
    let dq = tag[cstring](r, DEL_QUAL_TAG)
    if dq.isSome:
      return decodeQual(dq.get[i])
    else:
      return high(int)


proc delQualities(r: Record, q: var seq[uint8]): seq[uint8] =
  if len(q) != r.b.core.l_qseq:
    q.set_len(r.b.core.l_qseq)
  let dq = tag[cstring](r, DEL_QUAL_TAG)
  for i in 0..<int(r.b.core.l_qseq):
    if dq.isSome:
      q[i] = decodeQual(dq.get[i])
    else:
      q[i] = high(uint8) 


proc baseAlnQualityAt(r: Record, i: int): Natural =
  let baq = tag[cstring](r, BASE_ALN_QUAL_TAG)
  if baq.isSome:
    return decodeQual(baq.get[i])
  else:
    return high(int)


proc baseAlnQualities(r: Record, q: var seq[uint8]): seq[uint8] =
  if len(q) != r.b.core.l_qseq:
    q.set_len(r.b.core.l_qseq)
  let baq = tag[cstring](r, BASE_ALN_QUAL_TAG)
  for i in 0..<int(r.b.core.l_qseq):
    if baq.isSome:
      q[i] = decodeQual(baq.get[i])
    else:
      q[i] = high(uint8) 


proc insAlnQualityAt(r: Record, i: int): Natural =
  let iaq = tag[cstring](r, INS_ALN_QUAL_TAG)
  if iaq.isSome:
    return decodeQual(iaq.get[i])
  else:
    return high(uint8)


proc insAlnQualities(r: Record, q: var seq[uint8]): seq[uint8] =
  if len(q) != r.b.core.l_qseq:
    q.set_len(r.b.core.l_qseq)
  let iaq = tag[cstring](r, INS_ALN_QUAL_TAG)
  for i in 0..<int(r.b.core.l_qseq):
    if iaq.isSome:
      q[i] = decodeQual(iaq.get[i])
    else:
      q[i] = high(uint8)


proc delAlnQualityAt(r: Record, i: int): Natural =
  let daq = tag[cstring](r, DEL_ALN_QUAL_TAG)
  if daq.isSome:
    return decodeQual(daq.get[i])
  else:
    return high(uint8)


proc delAlnQualities(r: Record, q: var seq[uint8]): seq[uint8] =
  if len(q) != r.b.core.l_qseq:
    q.set_len(r.b.core.l_qseq)
  let daq = tag[cstring](r, DEL_ALN_QUAL_TAG)
  for i in 0..<int(r.b.core.l_qseq):
    if daq.isSome:
      q[i] = decodeQual(daq.get[i])
    else:
      q[i] = high(uint8)


proc getReadQualityBuffer*(r: Record): TReadQualityBuffer =
    result.baseQuals = newSeq[uint8]()# base qualities; already decoded
    result.baseAlnQuals = newSeq[uint8]()# base alignment qualities
    result.insQuals = newSeq[uint8]()# insertion qualities
    result.insAlnQuals = newSeq[uint8]()# insertion alignment qualities
    result.delQuals = newSeq[uint8]()# deletion qualities
    result.delAlnQuals = newSeq[uint8]()# deletion alignment qualities

    discard r.baseQualities(result.baseQuals)# from htsnim. see https://github.com/brentp/hts-nim/blob/13ccd2838a3d5e58991ff4cd50fe7ed41b434583/src/hts/bam.nim#L131
    discard r.baseAlnQualities(result.baseAlnQuals)
    discard r.insQualities(result.insQuals)
    discard r.insAlnQualities(result.insAlnQuals)
    discard r.delQualities(result.delQuals)
    discard r.delAlnQualities(result.delAlnQuals)



# FIXME use qual buffer instead
proc matchQual(r: Record, i: int, useMQ: bool): Natural =
  var q_m = high(int)
  let q_a = r.baseAlnQualityAt(i)
  let q_b = int(r.baseQualityAt(i))
  if useMQ:
    q_m = int(r.mapping_quality)
  mergeQuals(q_m, q_a, q_b)


# FIXME use qual buffer instead
proc insQual(r: Record, i: int, useMQ: bool): Natural =
  var q_m = high(int)
  let q_a = r.insAlnQualityAt(i)
  let q_i = r.insQualityAt(i)
  if useMQ:
    q_m = int(r.mapping_quality)
  mergeQuals(q_m, q_a, q_i)


# FIXME use qual buffer instead
proc delQual(r: Record, i: int, useMQ: bool): Natural =
  var q_m = high(int)
  let q_a = r.delAlnQualityAt(i)
  let q_d = r.delQualityAt(i)
  if useMQ:
    q_m = int(r.mapping_quality)
  mergeQuals(q_m, q_a, q_d)


proc newProcessor*[TStorage](storage: TStorage, useMQ: bool, minBQ: int):
  # minBQ = minimum base quality. everything below will be recorded as -1.
  # is later ignored/filtered by call(). 3 is default,
  # so that Illumina's Read Segment Quality Control Indicator" (#) gets ignored

  Processor[TStorage] {.inline.} =
  ## The default constructor for the 'Processor' type.
  Processor[TStorage](storage: storage,
                      matchQualityAt: matchQual,
                      insertionQualityAt: insQual,
                      deletionQualityAt: delQual,
                      minBQ: minBQ,
                      useMQ: useMQ)


proc processMatches*[TSequence](self: Processor,
                   readStart: int, refStart: int64, length: int,
                   read: Record, reference: TSequence,
                   nextevent: CigarElement) : void {.inline.} =
  ## Processes a matching substring between the read and the reference. All
  ## necessary information is available through the arguments. A matching
  ## substring consists of multiple contiguous matching bases.
  for offset in countUp(0, length - 1):
    let refOff = int(refStart + offset)# FIXME stupid
    let readOff = readStart + offset
    let bq = int(read.baseQualityAt(readOff))
    if bq >= self.minBQ:
      self.storage.recordMatch(refOff, $read.baseAt(readOff),
                                self.matchQualityAt(read, readOff, self.useMQ),
                                read.flag.reverse,
                                reference.baseAt(refOff))
    else:
      self.storage.recordMatch(refOff, $read.baseAt(readOff),
                                -1,# flag for later filtering
                                read.flag.reverse,
                                reference.baseAt(refOff))
    # Here we also need to record the indel qualities emitted from matches.
    # Just be careful to not count twice (hence check next op if at the end)
    if offset < length-1:
        self.storage.recordInsertion(refOff, $refSymbolAtIndel(read.flag.reverse),
          self.insertionQualityAt(read, readOff, self.useMQ),
          read.flag.reverse)
        self.storage.recordDeletion(refOff, $refSymbolAtIndel(read.flag.reverse),
          self.deletionQualityAt(read, readOff, self.useMQ),
          read.flag.reverse)
    elif nextevent.op == CigarOp.insert:
      self.storage.recordDeletion(refOff, $refSymbolAtIndel(read.flag.reverse),
        self.deletionQualityAt(read, readOff, self.useMQ),
        read.flag.reverse)
    elif nextevent.op == CigarOp.deletion:
      self.storage.recordInsertion(refOff, $refSymbolAtIndel(read.flag.reverse),
        self.insertionQualityAt(read, readOff, self.useMQ),
        read.flag.reverse)


proc processInsertion*[TSequence](self: Processor,
                     readStart: int, refIndex: int64, length: int,
                     read: Record, reference: TSequence): void {.inline.} =
  ## Processes an insertion on the read (wrt. the reference). All necessary
  ## information is available through the arguments. An insertion consists of
  ## one or more bases found on the read, but not on the reference.
  var value = ""
  for offset in countUp(readStart, readStart + length - 1):
    value &= read.baseAt(offset)

  # insertion is reported on the base that preceeds it
  self.storage.recordInsertion(refIndex - 1, value,
                               self.insertionQualityAt(read, readStart, self.useMQ),
                               read.flag.reverse)


proc processDeletion*[TSequence](self: Processor,
                    readIndex: int, refStart: int64, length: int,
                    read: Record, reference: TSequence): void {.inline.} =
  ## Processes an deletion on the read (wrt. the reference). All necessary
  ## information is available through the arguments. A deletion consists of one
  ## or more bases found on the read, but not on the reference.
  var value = ""
  for offset in countUp(refStart, refStart + length - 1):
    value &= reference.baseAt(int(offset))# stupid conversion
    self.storage.recordMatch(offset, $DEFAULT_BLANK_SYMBOL,
                             DEFAULT_BLANK_QUALITY,
                             read.flag.reverse,
                             reference.baseAt(int(offset)))# FIXME stupid int conversion

  # deletion is reported on the base that preceeds it
  self.storage.recordDeletion(refStart - 1, value,
                              self.deletionQualityAt(read, readIndex, self.useMQ),
                              read.flag.reverse)


proc beginRead*(self: Processor, read: Record): void {.inline.} =
  ## flush the storage up to the starting position.
  discard self.storage.flushUpTo(read.start)
  # buffer all read qualities for optimization (only parse qualities once)
  self.readQualityBuffer = read.getReadQualityBuffer()
     

proc done*(self: Processor): void {.inline.} =
  ## Finishes the processing, flushes the entire storage.
  discard self.storage.flushAll()
