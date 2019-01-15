# [EOF-Library](https://EOF-Library.com)
Libraries for coupling [Elmer FEM](https://www.csc.fi/web/elmer) and [OpenFOAM](https://openfoam.org/) + test cases. For detailed information how this software works check our [article](https://doi.org/10.1016/j.softx.2019.01.007):

Juris Vencels, Peter Råback, Vadims Geža,
__EOF-Library: Open-source Elmer FEM and OpenFOAM coupler for electromagnetics and fluid dynamics__,
SoftwareX, Volume 9, 2019, Pages 68-72, ISSN 2352-7110,
https://doi.org/10.1016/j.softx.2019.01.007.

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
docker run --rm -it -e HOST_USER_ID=$(id -u) -e HOST_USER_GID=$(id -g) -v ${PWD}:/home/openfoam/EOF-Library/runs eoflibrary/eof_elmer84_of6:latest
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
* Run simulation on 2 physical cores:
```
mpirun -n 2 mhdInterFoam -parallel : -n 2 ElmerSolver_mpi case.sif
```
* Simulation results will appear in host system *runs* folder


#### 2. Manual installation ####
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
* Update environment
```
. EOF-Library/etc/bashrc
```
* Check the MPI implementation and version (it is **important** that Elmer and OpenFOAM was compiled with the same version!)
```
which mpirun
mpirun --version
```
* Compile EOF-Library

```
eofCompile
```
* Colmpile OpenFOAM solver
```
cd EOF-Library/solvers/mhdInterFoam
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
