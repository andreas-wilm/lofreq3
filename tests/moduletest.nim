# standard library
import unittest
import os
import strutils
import osproc
# project specific
#/
# third party
#/


suite "module internal tests":
    let cmd = "nim c --stackTrace:on --lineTrace:on --checks:on --assertions:on --opt:none -r"
    for file in walkDirRec "../src/lofreqpkg/":
      if file.endswith(".nim"):
        test file:
          let (outp, errC) = execCmdEx(cmd & " " & file)
          check errC == 0
          var bin = file[0..^5]
          removeFile(bin)
