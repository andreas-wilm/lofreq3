## Implementation of an alternative procedure to samtools' pileup that does not
## require keeping pointers to all the reads in memory and provides support for
## streaming results to microservices as soon as possible. This was inspired by
## https://brentp.github.io/post/no-pile/.
##
## - Author: Filip Sodić <filip.sodic@gmail.com>
## - License: The MIT License

# standard
import times
import logging
import strutils
import sequtils

# third party
import hts

# project specific
import storage/slidingDeque
import recordFilter
import algorithm
import postprocessing
import ../region
import ../vcf
import ../call


var logger = newConsoleLogger(fmtStr = verboseFmtStr, useStderr = true)


proc fullPileup*(bamFname: string, faFname = "", regionsStr = "", bedFile = "",
                  handler: DataToVoid) : void =
  ## Performs the pileup over all chromosomes listed in the bam file.
  var bam: Bam
  var fai: Fai
  let numHTSReaderThreads = 1# see no improvement with 2 threads. likely all time spend on processing rather than unpacking
  if not open(bam, bamFname, index=true, threads=numHTSReaderThreads):
    quit("Could not open BAM file " & bamFname)

  if len(faFname)!=0:
    logger.log(lvlInfo, "Opening index for " & faFname)
    if not open(fai, faFname):
      quit("Could not open reference " & faFname)
  else:
    logger.log(lvlInfo, "No reference file given")


    for i in 0..<len(fai):
      let n = fai[i]
      let l = fai.chrom_len(fai[i])
      echo n & " " & $l
    #quit("FIXME")

  if len(regionsStr) != 0 and len(bedFile) != 0:
    quit("Can't read regions from bed and string at the same time")   

  var regions: seq[Region]
  if len(regionsStr) != 0:
    regions = toSeq(parseRegionsStr(regionsStr))
  elif len(bedFile) != 0:
    regions = toSeq(getBedRegions(bedFile))
  else:
    regions = toSeq(getBamRegions(bam))

  for reg in regions:
    logger.log(lvlInfo, "Starting pileup for " & $reg)

    var records = newRecordFilter(bam, reg.sq, reg.s, reg.e)

    let time = cpuTime()
    algorithm.pileup(fai, records, reg, handler)
    logger.log(lvlInfo, "Time taken to pileup reference ",
      reg.sq, " ", cpuTime() - time)


## "main" function. actually a pileup function with different postprocessing options
proc call*(bamFname: string, faFname: string, regions = "", bedFname = "",
           minVarQual: int = DEFAULT_MIN_VAR_QUAL,
           minAF: float = DEFAULT_MIN_AF,
           minCov: int = DEFAULT_MIN_COV,
           maxCov: int = DEFAULT_MAX_COV,
           minBQ: int = DEFAULT_MIN_BQ,
           noMQ: bool = not DEFAULT_USE_MQ,
           loglevel = 0, pileup = false, pretty = false) =

  if logLevel >= 3:
    setLogFilter(lvlDebug)
  elif logLevel == 2:
    setLogFilter(lvlInfo)
  elif logLevel == 1:
    setLogFilter(lvlNotice)
  elif logLevel == 0:
    setLogFilter(lvlWarn)
  else:
    quit("Invalid log level")

  var p: DataToVoid
  if pileup:
    if pretty:
      logger.log(lvlWarn, "Pretty printing is good for debugging,",
                 "but cannot be used for calling")
      p = toJsonAndPrettyPrint
    else:
      p = toJsonAndPrint
  else:
    if pretty:
      quit("Pretty print can only be used in conjuction with json")
    echo vcfHeader()
    callParams.minVarQual = minVarQual
    callParams.minAF  = minAF
    p = callAndPrint

  plpParams.minCov = minCov
  plpParams.maxCov = maxCov
  plpParams.minBQ = minBQ
  plpParams.useMQ = not noMQ
  fullPileup(bamFname, faFname, regions, bedFname, p)


