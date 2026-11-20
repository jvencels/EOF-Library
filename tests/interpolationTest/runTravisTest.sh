#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

nFails=0

elmerf90 -o InterpolationTester.so InterpolationTester.F90
wmake $EOF_HOME/solvers/interpolationTestFoam

ElmerGrid 2 2 meshElmer -metis 2 -nooverwrite
decomposePar

for case in *.sif; do
  rm TEST.PASSED*
  rm *.out
  if [ "0$FOAM_API" -ge 2312 ] 2>/dev/null && [ -z "$EOF_COMMSPLIT" ]; then
    # openfoam.com with the installed libPstream: ElmerSolver first, OpenFOAM split by application number
    mpirun -n 2 ElmerSolver_mpi $case : -n 2 interpolationTestFoam -parallel -mpi-split-by-appnum
  else
    mpirun -n 2 interpolationTestFoam -parallel : -n 2 ElmerSolver_mpi $case
  fi
  if [[ $(< TEST.PASSED*) != "1" ]]; then
    nFails=$((nFails+1))
  fi
done

exit $nFails
