#!/bin/bash
cd "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

nFails=0

wmake $EOF_HOME/solvers/mhdInterFoamv2606

ElmerGrid 2 2 meshElmer -metis 2
setFields
decomposePar

for case in *.sif; do
  rm TEST.PASSED*
  if [ -z "$EOF_COMMSPLIT" ]; then
    # Installed libPstream: ElmerSolver first, OpenFOAM split by application number
    mpirun -n 2 ElmerSolver_mpi $case : -n 2 mhdInterFoamv2606 -parallel -mpi-split-by-appnum
  else
    # commSplit libPstream override: any order
    mpirun -n 2 mhdInterFoamv2606 -parallel : -n 2 ElmerSolver_mpi $case
  fi
  if [[ $(< TEST.PASSED*) != "1" ]]; then
    nFails=$((nFails+1))
  fi
done

exit $nFails
