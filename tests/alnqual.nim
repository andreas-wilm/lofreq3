# standard library
import unittest
import math
import osproc
#import strutils
import tempfile
import strformat

# project specific
import ../src/lofreqpkg/alnqual
#import ../src/lofreqpkg/utils

# third party
import hts


let lofreq = "../lofreq"


suite "alignment qualities (baq, ai, ad)":

  test "simple":
    const fa = "alnqual/ref.fa"
    const bam = "alnqual/alignments.bam"

    # make sure call works and  nothing is predicted
    var tmpfd: File
    var tmpname: string
    (tmpfd, tmpname) = mkstemp()
    tmpfd.close
    let cmd = fmt"{lofreq} alnqual -f {fa} -b {bam} > {tmpname}"
    echo "Testing: " & cmd
    let outp = execProcess(cmd)

    var b:Bam
    open(b, tmpname, index=false, fai=fa)
    for record in b:# only one record in there
      let ai = tag[string](record, AI_TAG)
      check not ai.isNone
      check ai.get == "~~~!~~~"

      let ad = tag[string](record, AD_TAG)
      check not ad.isNone
      check ad.get == "!~~~~~~"
      
      let baq = tag[string](record, BAQ_TAG)
      check not baq.isNone
      check baq.get == "!!!!55!"
      
      #  ../lofreq.git/src/lofreq/lofreq alnqual $bam $fa
      #@HD     VN:1.4  SO:coordinate
      #@SQ     SN:ref  LN:14
    # 1       0       ref     1       60      1M2D3M2I1M      *       0       0       AAAAGCC 5555555 BI:Z:???????    BD:Z:IIIIIII lb:Z:!!!!55!     ai:Z:~~~!~~~    ad:Z:!~~~~~~

