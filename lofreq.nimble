version       = "3.0.0"
author        = "Andreas Wilm and others"
description   = "LoFreq Reimplementation in Nim"
license       = "MIT"

requires "nim >= 1.0", "cligen >= 0.9.16", "hts >= 0.3.3", "tempfile >= 0.1.5"
# know to work with
#  nim = 0.19.4, cligen = 0.9.19, hts = 0.2.8, tempfile = 0.1.7


srcDir = "src"

bin = @["lofreq", "vcfeval"]

skipDirs = @["tests"]
skipExt = @["nim"]

task test, "run tests":
  withDir "tests":
    exec "nim c --lineDir:on --debuginfo -r all"

