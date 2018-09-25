version       = "3.0.0"
author        = "Andreas Wilm and others"
description   = "LoFreq Reimplementation in Nim"
license       = "MIT"

requires "nim >= 0.18", "cligen >= 0.9.11", "hts >= 0.2.4"

srcDir = "src"

bin = @["lofreq"]

skipDirs = @["tests"]
skipExt = @["nim"]

task test, "run tests":
  withDir "tests":
    exec "nim c --lineDir:on --debuginfo -r all"

