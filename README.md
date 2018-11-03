# [EOF-Library](https://EOF-Library.com) // [![Build Status](https://travis-ci.org/jvencels/EOF-Library.svg?branch=master)](https://travis-ci.org/jvencels/EOF-Library)
Libraries for coupling [Elmer FEM](https://www.csc.fi/web/elmer) and [OpenFOAM](https://openfoam.org/) + test cases. For detailed information how this software works check our [article preprint](http://dx.doi.org/10.13140/RG.2.2.12907.39203).

## About ##
This software is maintained by [Juris Vencels](https://lv.linkedin.com/in/vencels) from *University of Latvia* in cooperation with *CSC - IT Center for Science Ltd.* (Finland).

* Software improvements and suggestions are very welcome. 
* If you find this software useful then drop me a message, we could promote each other work.
* The best way to support development of this library is buying consulting services.
* GPL v3 license allows commercial use of this software. 

## Introduction ##
We use EOF-Library for solving magnetohydrodynamic (MHD) problems, e.g., electromagnetic induction melting, stirring, levitation and other cases where metals interact with electromagnetic fields. 

Other users have found it useful for plasma physics, convective cooling of electrical devices and other applications. EOF-Library allows coupling internal fields between almost any Elmer and OpenFOAM solvers.

## Requirements ##
* **Both Elmer and OpenFOAM must use the same OpenMPI version!**

## How to ##
There are two options to install and use this software:
1. __Docker__ install (best for *beginners* and running on *clouds*) - **Linux, Windows, MacOS**
2. __Complete__ install (best for *developers* and *advanced* users) - **Linux**

Performance-wise these two options are comparable. Docker installation comes with OpenFOAM and Elmer installed, and environment is set. This is preferred option for most users and for running simulations on multiple computers or cloud.

On the other hand complete installation gives users more flexibility, it is preferred option for developers who want to work on their own solvers or have full control over software and its source code.

#### 1. Docker ####
First, you will need to install docker on your OS:
1. **Ubuntu/Linux** (preferred) - https://docs.docker.com/install/linux/docker-ce/ubuntu/
2. **Windows** - https://docs.docker.com/docker-for-windows/install/
https://docs.docker.com/docker-for-windows/install/
3. **MacOS** - https://docs.docker.com/docker-for-mac/install/


Then, follow commands below to install the software & run demo simulation.
* Create an empty folder
```
mkdir runs
cd runs
```
* Run Docker image and bind mount current host system folder *${PWD}* to newly created *EOF-Library/runs* folder
```
docker run -it -v ${PWD}:/home/openfoam/EOF-Library/runs eoflibrary/eof_of6:latest
```
* Update EOF-Library and compile it
```
eofUpdate
```
* Compile OpenFOAM solver
```
cd EOF-Library
wmake solvers/mhdInterFoam6
```
* Copy test simulation
``` 
cp -r tests/levitation2D runs
```
* Prepare case
```
cd runs/levitation2D
setFields
decomposePar
ElmerGrid 2 2 meshElmer -metis 2
```
* Run simulaiton on 2 physical cores
```
mpirun --allow-run-as-root -n 2 mhdInterFoam -parallel : -n 2 ElmerSolver_mpi case.sif
```
* Simulation results will appear in host system *runs* folder


#### 2. Complete installation ####
* Get `git`, `cmake`, `gfortran`, `blas` and `lapack`.
```
sudo apt-get install git cmake gfortran libblas-dev liblapack-dev
```
* Download OpenFOAM 6.0 source code, configure and compile it (https://openfoam.org/download/6-source/)
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
* Compile libraries

```
cd EOF-Library/libs
elmerf90 -o Elmer2OpenFOAM.so Elmer2OpenFOAM.F90
elmerf90 -o OpenFOAM2Elmer.so OpenFOAM2Elmer.F90
```
* Add this line to `.bashrc`
```
export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH
```
* Colmpile OpenFOAM solver
```
cd coupleElmer
wmake
cd ../../solvers/mhdInterFoam
wmake
```

* Copy test

```
cd ../..
mkdir runs
cd runs
cp -r ../tests/levitation2D .
```

* Prepare test

```
cd levitation2D
setFields
decomposePar
ElmerGrid 2 2 meshElmer -metis 2
```

* Run OpenFOAM on 2 processes and Elmer on 2 processes

```
mpirun -np 2 mhdInterFoam -parallel : -np 2 ElmerSolver_mpi
```
