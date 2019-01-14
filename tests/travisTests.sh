#!/bin/bash

cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. /opt/openfoam6/etc/bashrc
. ../etc/bashrc

eofCompile &>> travisTests.log

nTestsPass=0
nTestsFail=0
 
for testDir in */ ; do
  echo "Running test:" $testDir
  ./$testDir/runTravisTest.sh &> tmp.log
  status=$?
  if [[ $status -eq 0 ]] ; then
    echo "  Pass!"
    nTestsPass=$((nTestsPass+1))
  else
    echo "  Fail! $status"
    cat tmp.log >> travisTests.log
    nTestsFail=$((nTestsFail+status))
  fi
done

if [ $nTestsFail -gt 0 ] ; then
  echo "$nTestsFail tests failed!"
  echo "==================================="
  cat travisTests.log
  echo "==================================="
  echo "$nTestsFail tests failed!"
  exit 1
else
  echo "All tests passed!"
  exit 0
fi
