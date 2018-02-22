# EOF-Library (https://EOF-Library.com)
Libraries for coupling Elmer and OpenFOAM + test cases. Currently used for solving coupled problems, e.g., electromagnetic induction heating and MHD with free surface, but users can use this library for coupling any Elmer and OpenFOAM solvers.

This software is maintained by the *Laboratory for mathematical modelling of environmental and technological processes* (University of Latvia) in cooperation with *CSC - IT Center for Science Ltd.* (Finland).

## Requirements ##
* Tested on Ubuntu 14.04 and Ubuntu 16.04.
* Currently supported **OpenFOAM** version **5.0** compiled from source (https://openfoam.org/download/5-0-source/). Other OpenFOAM versions currently are not supported.
* Both Elmer and OpenFOAM must use the same OpenMPI version.
* You will need `git`, `cmake`, `gfortran`, `blas` and `lapack`. 

```
sudo apt-get install git cmake gfortran libblas-dev liblapack-dev
```

## How to ##

* Download OpenFOAM 5.0 source code, configure and compile it (https://openfoam.org/download/5-0-source/)
* Download Elmer from developers repo (https://github.com/ElmerCSC/elmerfem).

* Configure and compile Elmer with `-DWITH_MPI=TRUE` by following these steps (https://www.csc.fi/web/elmer/sources-and-compilation). In short .., we create build folder and call CMake from this folder:
```
cmake -DWITH_MPI=TRUE -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../mpi-install ../elmerfem
make -j install
```

* Download EOF-Library:

```
git clone https://github.com/jvencels/EOF-Library
```
* Add this line to `.bashrc`
```
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
```
* Check the MPI implementation and version (it is **important** that Elmer and OpenFOAM was compiled with the same version!)
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
ElmerGrid 2 2 mesh_Elmer -metis 2
```

* Copy `Elmer2OpenFOAM.so`, `Elmer2OpenFOAM.mod`, `OpenFOAM2Elmer.so` and `OpenFOAM2Elmer.mod` to this folder

* Run OpenFOAM on 2 processes and Elmer on 2 processes

```
mpirun -np 2 mhdInterFoam -parallel : -np 2 ElmerSolver_mpi
```

* Postprocessing

```
reconstructPar
paraFoam
```
