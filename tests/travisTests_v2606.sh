#!/bin/bash
# Run the tests that support openfoam.com v2606. Set OPENFOAM_BASHRC to the etc/bashrc
# of the OpenFOAM installation (default: the Debian/Ubuntu package location). Without
# EOF_COMMSPLIT the installed libPstream is used and ElmerSolver must start first on the
# mpirun line with OpenFOAM given -mpi-split-by-appnum; with EOF_COMMSPLIT=1 the commSplit
# libPstream override is built and either order works.

cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. ${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2606/etc/bashrc}
. ../etc/bashrc

eofCompile &>> travisTests_v2606.log

tests="interpolationTest levitation2Dv2606"

nTestsPass=0
nTestsFail=0

for testDir in $tests ; do
  echo "Running test:" $testDir
  ./$testDir/runTravisTest.sh &> tmp.log
  status=$?
  if [[ $status -eq 0 ]] ; then
    echo "  Pass!"
    nTestsPass=$((nTestsPass+1))
  else
    echo "  Fail! $status"
    tail -500 tmp.log >> travisTests_v2606.log
    nTestsFail=$((nTestsFail+status))
  fi
done

if [ $nTestsFail -gt 0 ] ; then
  echo "$nTestsFail tests failed!"
  echo "==================================="
  cat travisTests_v2606.log
  echo "==================================="
  echo "$nTestsFail tests failed!"
  exit 1
else
  echo "All tests passed!"
  exit 0
fi
