# standard library
import unittest
import osproc
import tempfile

# project specific
#/

# third party
import hts/vcf


let lofreq = "../lofreq"


suite "coverage":

  test "Checking mincov":
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    let cmd = lofreq & " call -f call_samples/cov.fa -b call_samples/cov.bam --mincov 4 > " & tmpname
    #echo "Testing: " & cmd
    discard execProcess(cmd)
    var v:VCF
    discard open(v, tmpname)
    var nvars = 0
    for rec in v:
        inc nvars
    check nvars == 0


  test "Checking maxcov":
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    let cmd = lofreq & " call -f call_samples/cov.fa -b call_samples/cov.bam --maxcov 1 > " & tmpname
    #echo "Testing: " & cmd
    discard execProcess(cmd)
    var v:VCF
    discard open(v, tmpname)
    var nvars = 0
    for rec in v:
        inc nvars
    check nvars == 2


  # simple test, mainly an excuse to test htsnim vcf
  test "DP and DP4":
    # from lofreq 2.1.5
    let DP_REF = @[1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1]
    let DP4_REF = @[
      @[0,0,1,0],
      @[0,0,2,0],
      @[0,0,3,0],
      @[2,0,1,0],
      @[0,0,3,0],
      @[0,0,3,0],
      @[0,0,3,0],
      @[0,0,3,0],
      @[2,0,1,0],
      @[0,0,2,0],
      @[0,0,3,0],
      @[0,0,1,0],
      ]
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    let cmd = lofreq & " call -f call_samples/cov.fa -b call_samples/cov.bam > " & tmpname
    #echo "Testing: " & cmd
    let outp = execProcess(cmd)
    var v:VCF
    discard open(v, tmpname)
    var nvars = 0
    for rec in v:
      var dps = newSeq[int32](1)
      var dp4 = newSeq[int32](4)
      check rec.info.get("DP", dps) == Status.OK
      check dps[0] == DP_REF[nvars]
      check rec.info.get("DP4", dp4) == Status.OK
      var dp4_ref = newSeq[int32](4)
      for i in 0..3:
        dp4_ref[i] = int32(DP4_REF[nvars][i])
      #echo "FIXME : " & $rec & " ref " & $dp4_ref & " is " & $dp4
      check dp4 == dp4_ref
      inc nvars


