# OpenFOAM_Elmer
Libraries for OpenFOAM and Elmer MPI coupling + test cases. Used for solving coupled electromagnetic induction and MHD problems. This software is developed at Laboratory for mathematical modelling of environmental and technological processes (University of Latvia). 

## Requirements ##

Tested on:
* Ubuntu 14.04
* Open MPI v2.0.1 (https://www.open-mpi.org/software/ompi/v2.0/)

## How to ##
* Download Elmer (https://github.com/ElmerCSC/elmerfem.git)
* Compile Elmer with `-DWITH_MPI=TRUE`
* Compile `Elmer2OpenFOAM.F90` and `OpenFOAM2Elmer.F90`

```
elmerf90 -o Elmer2OpenFOAM.so Elmer2OpenFOAM.F90
elmerf90 -o OpenFOAM2Elmer.so OpenFOAM2Elmer.F90
```

* Download OpenFOAM (http://openfoam.org/) source and compile
* Download modified OpenFOAM libraries, solvers and tests:

```
git clone https://github.com/jvencels/OpenFOAM_Elmer.git
```

* Compile library and solver

```
cd OpenFOAM_Elmer/libs/commSplit
wmake
cd ../coupleElmer
wmake
cd ../mhdInterFoam
wmake
```
* Set environment variable (add to `.bashrc`)

```
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
```

* Copy test

```
cd ../..
mkdir runs
cd runs
cp -r ../tests/mhdInterFoamTest_2D .
```

* Prepare test

```
cd mhdInterFoamTest_2D
decomposePar
```

* Run OpenFOAM on 2 processes and Elmer on 1 process

```
mpirun -np 2 mhdInterFoam -parallel : -np 1 ElmerSolver_mpi
```

* Postprocessing

```
reconstructPar
paraFoam
```

* To open Elmer results, open `mesh_Elmer/case_0001.vtu`
