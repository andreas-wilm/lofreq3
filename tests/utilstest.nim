import unittest
import math


import ../src/lofreqpkg/utils

    
suite "qual and prob conversion":
  #setup:
  #   echo "suite setup: run before each test"

  #teardown:
  #  echo "suite teardown: run after each test"

  test "prob2qual":
    check prob2qual(0.05) == 13
    check prob2qual(0.01) == 20
    check prob2qual(0.001) == 30

  test "qual2prob":
    check abs(qual2prob(13) - 0.05) < 0.005
    check abs(qual2prob(20) - 0.01) < 0.001
    check abs(qual2prob(30) - 0.001) < 0.0001


