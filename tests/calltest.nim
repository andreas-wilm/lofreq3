# standard library
import unittest
import osproc
import tempfile

# project specific
import ../src/lofreqpkg/call

# third party
import hts/vcf


let lofreq = "../lofreq"


suite "call":


  # issue 14: fails to remove zombie tests before calling
  test "remove zombie entries from plp before calling":
    # make sure call works and  nothing is predicted
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    let cmd = lofreq & " call_from_plp -p call_samples/zombie.json > " & tmpname
    #echo "Testing: " & cmd
    discard execProcess(cmd)
    var v:VCF
    discard open(v, tmpname)
    var nvars = 0
    for rec in v:
        inc nvars
        
    check nvars == 0


  # simple test, mainly an excuse to test htsnim vcf
  test "simple-vars":
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    #let plp_cmd = lofreq & " pileup -b call_samples/simple-vars.bam -f call_samples/NC_000913.n200.fa -r NC_000913:1-200 --noMQ"
    #let call_cmd = lofreq & " call -p - -m 13"
    #let cmd = plp_cmd & " | " & call_cmd & " > " & tmpname
    let cmd = lofreq & " call -b call_samples/simple-vars.bam -v 13 -f call_samples/NC_000913.n200.fa -r NC_000913:1-200 --noMQ > " & tmpname
    #echo "Testing: " & cmd
    let outp = execProcess(cmd)
    var v:VCF
    discard open(v, tmpname)
    var nvars = 0
    for rec in v:
        inc nvars
        var dps = new_seq[int32](1)
        var sbs = new_seq[int32](1)
        var afs = new_seq[float32](1)
        var dp4 = new_seq[int32](4)
        check rec.info.get("DP", dps) == Status.OK
        check dps[0] == 2
        check rec.info.get("SB", sbs) == Status.OK
        check sbs[0] == 0
        check rec.info.get("AF", afs) == Status.OK
        check abs(afs[0] - 0.5) < 0.0000001
        check rec.info.get("DP4", dp4) == Status.OK
        # don't know how to check DP4
        #check dp4 == @[0'i32, 1'i32, 1'i32, 0'i32]
        check rec.QUAL == 20

    check nvars == 10

  test "direct call vs call from plp":
    var
      tmpfd1: File
      tmpfd2: File
    var
      tmpname1: string
      tmpname2: string
    var
      output: TaintedString
      exitCode: int

    (tmpfd1, tmpname1) = mkstemp()
    tmpfd1.close
    (tmpfd2, tmpname2) = mkstemp()
    tmpfd2.close

    let cmd1 = lofreq & " call -b call_samples/simple-vars.bam -f call_samples/NC_000913.n200.fa -r NC_000913:1-200 > " & tmpname1
    #echo "Testing: " & cmd1
    (output, exitCode) = execCmdEx(cmd1)
    check(exitCode == 0)

    let plp_cmd = lofreq & " call -b call_samples/simple-vars.bam -f call_samples/NC_000913.n200.fa -r NC_000913:1-200 -p"
    let call_cmd = lofreq & " call_from_plp -p -"
    let cmd2 = plp_cmd & " | " & call_cmd & " > " & tmpname2
    #echo "Testing: " & cmd2
    (output, exitCode) = execCmdEx(cmd2)
    check(exitCode == 0)

    let diff_cmd = "diff -q " & tmpname1 & " " & tmpname2
    (output, exitCode) = execCmdEx(diff_cmd)
    #echo "Diff command: " & diff_cmd
    check(exitCode == 0)

