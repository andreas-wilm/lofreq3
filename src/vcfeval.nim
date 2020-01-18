## Simple bcftools isec based vcf evaluator
##
# standard
import cligen
import os
import math
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


proc vcfEvalMain(vcfTruth, vcfTest: string, varType = "") =
  var varTypes: seq[string]

  if len(varType) > 0:
    if varType in allowedVarTypes:
      varTypes.add(varType)
    else:
      quit("Invalid variant type {varType}. Must be one of {allowedVarTypes}")
  else:
    varTypes = allowedVarTypes

  for v in varTypes:
    let classes = vcfEval(vcfTruth, vcfTest, v)
    var recall = classes.tp / (classes.tp + classes.fn)
    if recall.classify == fcNaN:
      recall = -1.0
    var precision = classes.tp / (classes.tp + classes.fp)
    if precision.classify == fcNAN:
      precision = -1.0
    echo fmt"{v}: TP={classes.tp} FP={classes.fp} FN={classes.fn}"
    echo fmt"{v}: recall={recall:.2f}"
    echo fmt"{v}: precision={precision:.2f}"


when isMainModule:
  import cligen
  dispatch(vcfEvalMain,
    short = {"vcfTruth": 'r', "vcfTest": 't', "varType": 'v'},
    help = {"vcfTruth": "Truth VCF", "vcfTest": "Test VCF", "varType": "Variant type"})
