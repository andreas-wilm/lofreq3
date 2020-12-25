## Simple bcftools isec based vcf evaluator
##
# standard
import cligen
import os
#import math
import osproc
import strformat
# project specific
# #/
# third party
import tempfile
import hts


const allowedVarTypes = @["snp", "indel"]
# FIXME probably a good idea to normalize first


type classes = object
  tp: int
  fp: int
  fn: int


proc vcfNumVars(vcfFile: string): int =
  var vcfFH: VCF
  doAssert open(vcfFH, vcfFile)
  for variant in vcfFH:
    inc result


proc vcfEval(vcfTruth, vcfTest, varType: string): classes =
  doAssert varType in allowedVarTypes
  var odir = mkdtemp()
  let cmd = fmt"bcftools isec -p {odir} --include ""TYPE='{vartype}'"" {vcftruth} {vcftest}"
  let (output, errCode) = execCmdEx(cmd)
  if errCode != 0:
    quit(fmt"Command '{cmd}' failed with exit code {errCode} and following error message: {output}")

  result.tp = vcfNumVars(joinPath(odir, "0002.vcf"))
  result.fp = vcfNumVars(joinPath(odir, "0001.vcf"))
  result.fn = vcfNumVars(joinPath(odir, "0000.vcf"))

  echo fmt"{vartype}: Keeping bcftools isec tmpdir {odir}"


proc vcfEvalMain(vcfTruth, vcfTest: string, varType = "",
  minRecall = 0.0, minPrecision = 0.0, minTP = 0, maxFP = -1, maxFN = -1) =
  var varTypes: seq[string]

  if len(varType) > 0:
    if varType in allowedVarTypes:
      varTypes.add(varType)
    else:
      let allowedVarTypesStr = allowedVarTypes.join(" ")
      quit(fmt"Invalid variant type {varType}. Must be one of {allowedVarTypesStr}")
  else:
    varTypes = allowedVarTypes

  for v in varTypes:
    let classes = vcfEval(vcfTruth, vcfTest, v)
    var recall = classes.tp / (classes.tp + classes.fn)
    #if recall.classify == fcNaN:
    #  recall = -1.0
    var precision = classes.tp / (classes.tp + classes.fp)
    #if precision.classify == fcNAN:
    #  precision = -1.0

    # fail if not passing required tests
    if classes.tp < minTP:
      quit(fmt"{v}: TP={classes.tp}<{minTP}")
    if maxFP != -1 and classes.fp > maxFP:
       quit(fmt"{v}: FP={classes.fp}>{maxfp}")
    if maxFN != -1 and classes.fn > maxFn:
      quit(fmt"{v}: FP={classes.fn}>{maxfn}")
    if recall < minRecall:
      quit(fmt"{v}: recall={recall}<{minRecall}")
    if precision < minPrecision:
      quit(fmt"{v}: precision={precision}<{minPrecision}")

    # report values
    echo fmt"{v}: TP={classes.tp}"
    echo fmt"{v}: FP={classes.fp}"
    echo fmt"{v}: FN={classes.fn}"
    echo fmt"{v}: recall={recall:.2f}"
    echo fmt"{v}: precision={precision:.2f}"


when isMainModule:
  import cligen
  dispatch(vcfEvalMain,
    short = {"vcfTruth": 'r', "vcfTest": 't', "varType": 'v'},
    help = {"vcfTruth": "Truth VCF", "vcfTest": "Test VCF", "varType": "Variant type",
    "minRecall": "Fail if recall is below this value",
    "minPrecision": "Fail is precision is below this vaule",
    "minTP": "Fail is TP is below this value",
    "maxFP": "Fail if FP is higher than this value (neg. value == ignored)",
    "maxFN": "Fail is FN is higher than this value (neg. value == ignored)"})

