#!/bin/bash

source /opt/openfoam5/etc/bashrc
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH

ElmerGrid 2 2 meshElmer -metis 2 -nooverwrite &> logElmerGrid.out
decomposePar -force &> logDecomposePar.out

#mpirun --allow-run-as-root -n 2 interpolationTestFoam -parallel : -n 2 ElmerSolver_mpi "elem.sif"
mpirun --allow-run-as-root -n 2 interpolationTestFoam -parallel : -n 2 ElmerSolver_mpi "dg.sif"
#mpirun --allow-run-as-root -n 2 interpolationTestFoam -parallel : -n 2 ElmerSolver_mpi "ip.sif"
