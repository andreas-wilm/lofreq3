# standard library
import unittest
import osproc
import tempfile
import strutils

# project specific
import ../src/lofreqpkg/call

# third party
import hts/vcf


let lofreq = "../lofreq"


suite "bugs":


  # should not call variants on non-ACGT positions
  test "bug18-ignore-nonACGT-in-ref":
    # make sure call works and  nothing is predicted
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    let reffa="bugs/bug18-ignore-nonACGT-in-ref/ref.fa"
    let bam="bugs/bug18-ignore-nonACGT-in-ref/alignments.bam"
    let cmd = lofreq & " call -f " & reffa & " -b " & bam & " > " & tmpname
    #echo "Testing: " & cmd
    discard execProcess(cmd)
    var v:VCF
    discard open(v, tmpname)
    var nvars = 0
    for rec in v:
        inc nvars
        check rec.ALT[0] in "ACGT"    
    check nvars == 3# only calling on positions CGT



