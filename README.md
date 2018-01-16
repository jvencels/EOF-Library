# EOF-Library (http://EOF-Library.com)
Libraries for coupling Elmer and OpenFOAM + test cases. Used for solving coupled problems, e.g., electromagnetic induction and MHD. This software is maintained by the *Laboratory for mathematical modelling of environmental and technological processes* (University of Latvia) in cooperation with *CSC - IT Center for Science Ltd.* (Finland).

## Requirements ##
* Tested on Ubuntu 14.04 and Ubuntu 16.04.
* Currently supported **OpenFOAM** version **5.0**. We suggest to compile it from the source code since compiled version from repositories can be incompatible. 
* Both Elmer and OpenFOAM must use the same OpenMPI version
* You will need `git`, `cmake`, `gfortran`, `blas` and `lapack`. 

```
sudo apt-get install git cmake gfortran libblas-dev liblapack-dev
```

## How to ##

* Download OpenFOAM 5.0 source code, configure and compile it (https://openfoam.org/download/5-0-source/)
* Download Elmer (https://github.com/jvencels/elmerfem.git). You also can obtain it from developers repo (https://github.com/ElmerCSC/elmerfem) without any guarantee that coupling will work. Alternately, Elmer could be installed from launchpad (https://launchpad.net/~elmer-csc-ubuntu/+archive/ubuntu/elmer-csc-ppa). This has been verified to work in Ubuntu 16.04 with EOF-Library software.

* Configure and compile Elmer with `-DWITH_MPI=TRUE` by following these steps (https://www.csc.fi/web/elmer/sources-and-compilation).

* Download modified OpenFOAM libraries, solvers and tests:

```
git clone https://github.com/jvencels/EOF-Library
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
cd EOF-Library/libs/commSplit
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
