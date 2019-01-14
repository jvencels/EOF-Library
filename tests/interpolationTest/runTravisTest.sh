#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

nFails=0

elmerf90 -o InterpolationTester.so InterpolationTester.F90
wmake $EOF_HOME/solvers/interpolationTestFoam

ElmerGrid 2 2 meshElmer -metis 2 -nooverwrite
decomposePar

for case in *.sif; do
  rm TEST.PASSED*
  mpirun -n 2 interpolationTestFoam -parallel : -n 2 ElmerSolver_mpi $case
  if [[ $(< TEST.PASSED*) != "1" ]]; then
    nFails=$((nFails+1))
  fi
done

exit $nFails
