#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

nFails=0

wmake $EOF_HOME/solvers/mhdInterFoam6

ElmerGrid 2 2 meshElmer -metis 2
setFields
decomposePar

for case in *.sif; do
  rm TEST.PASSED*
  mpirun -n 2 mhdInterFoam -parallel : -n 2 ElmerSolver_mpi $case
  if [[ $(< TEST.PASSED*) != "1" ]]; then
    nFails=$((nFails+1))
  fi
done

exit $nFails
