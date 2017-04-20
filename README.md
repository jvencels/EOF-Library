# OpenFOAM_Elmer
Libraries for OpenFOAM and Elmer MPI coupling + test cases. Used for solving coupled electromagnetic induction and MHD problems. This software is developed at Laboratory for mathematical modelling of environmental and technological processes (University of Latvia). 

## Requirements ##
Tested on Ubuntu 14.04 and Ubuntu 16.04.
You will need `git`, `cmake`, `gfortran`, `blas` and `lapack`. 

```
sudo apt-get install git cmake gfortran libblas-dev liblapack-dev
```

## How to ##

* Download, install and configure OpenFOAM (http://openfoam.org/download/4-1-ubuntu/)
* Download Elmer (https://github.com/jvencels/elmerfem.git). You also can obtain it from developers repo (https://github.com/ElmerCSC/elmerfem) without any guarantee that coupling will work. Alternately, Elmer could be installed from launchpad (https://launchpad.net/~elmer-csc-ubuntu/+archive/ubuntu/elmer-csc-ppa). This has been verified to work in Ubuntu 16.04 with OpenFOAM_Elmer software.

* Configure and compile Elmer with `-DWITH_MPI=TRUE` by following these steps (https://www.csc.fi/web/elmer/sources-and-compilation).

* Download modified OpenFOAM libraries, solvers and tests:

```
git clone https://github.com/jvencels/OpenFOAM_Elmer.git
```
* Add this line to `.bashrc`
```
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
```
* Check the MPI implementation and version
```
which mpirun
mpirun --version
```
* Compile libraries and solver

```
cd OpenFOAM_Elmer/libs/commSplit
wmake
cd ../coupleElmer
wmake
cd ../mhdInterFoam
wmake
cd ../Elmer2OpenFOAM
elmerf90 -o Elmer2OpenFOAM.so Elmer2OpenFOAM.F90
cd ../OpenFOAM2Elmer
elmerf90 -o OpenFOAM2Elmer.so OpenFOAM2Elmer.F90
```

* Add this line to `.bashrc`


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
setFields
decomposePar
```

* Copy `Elmer2OpenFOAM.so`, `Elmer2OpenFOAM.mod`, `OpenFOAM2Elmer.so` and `OpenFOAM2Elmer.mod` to this folder

* Run OpenFOAM on 2 processes and Elmer on 1 process

```
mpirun -np 2 mhdInterFoam -parallel : -np 1 ElmerSolver_mpi
```

* Postprocessing

```
reconstructPar
paraFoam
```
