version       = "3.0.0"
author        = "Andreas Wilm and others"
description   = "LoFreq Reimplementation in Nim"
license       = "MIT"

requires "nim >= 0.19", "cligen >= 0.9.16", "hts >= 0.2.5", "tempfile >= 0.1.5"

srcDir = "src"

bin = @["lofreq"]

skipDirs = @["tests"]
skipExt = @["nim"]

task test, "run tests":
  withDir "tests":
    exec "nim c --lineDir:on --debuginfo -r all"

