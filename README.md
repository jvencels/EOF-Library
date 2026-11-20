# [EOF-Library](https://EOF-Library.com)
Libraries for coupling [Elmer FEM](https://www.csc.fi/web/elmer) and [OpenFOAM](https://openfoam.org/) + test cases. For detailed information how this software works check our [article](https://doi.org/10.1016/j.softx.2019.01.007).

___If you use EOF-Library, please cite:___
Juris Vencels, Peter Råback, Vadims Geža,
_EOF-Library: Open-source Elmer FEM and OpenFOAM coupler for electromagnetics and fluid dynamics_,
SoftwareX, Volume 9, 2019, Pages 68-72, ISSN 2352-7110,
https://doi.org/10.1016/j.softx.2019.01.007.

## About ##
This software is maintained by [Juris Vencels](https://lv.linkedin.com/in/vencels) from *EOF Consulting LLC* in cooperation with *CSC - IT Center for Science Ltd.* (Finland).

* Software improvements and suggestions are very welcome. 
* If you find this software useful then drop me a message.
* The best way to support development of this library is buying consulting services.
* GPL v3 license allows commercial use of this software.

## Introduction ##
EOF-Library couples internal fields between Elmer FEM and OpenFOAM. Applications are 
* Magnetohydrodynamics (MHD)
* Microwave ehating
* Plasma physics
* Convective cooling of electrical devices and machines

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
* Download OpenFOAM source code, configure and compile it
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
cd EOF-Library/solvers/mhdInterFoam6
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

#### 3. openfoam.com (v1812, and v2312 or newer) ####

OpenFOAM from openfoam.com is supported in two flavours. In both cases source
the OpenFOAM `etc/bashrc` and `EOF-Library/etc/bashrc`, put the Elmer `bin`
directory on the `PATH`, and run `eofCompile`; the solver for these versions is
`solvers/mhdInterFoamv2606` and the test case is `tests/levitation2Dv2606`.

* **v1812**: `eofCompile` builds the `libs/commSplit` override of `libPstream`
  first, as for openfoam.org versions below 6. The job is started in the usual
  order:
```
mpirun -np 2 mhdInterFoamv2606 -parallel : -np 2 ElmerSolver_mpi case.sif
```

* **v2312 or newer, default (native world split)**: the installed `libPstream`
  is used and OpenFOAM splits `MPI_COMM_WORLD` itself with its run-time option
  `-mpi-split-by-appnum`. ElmerSolver must then be the **first** application on
  the `mpirun` line:
```
mpirun -np 2 ElmerSolver_mpi case.sif : -np 2 mhdInterFoamv2606 -parallel -mpi-split-by-appnum
```
  Starting the job in the old order makes the coupler stop with the message
  "No partner ranks".

* **v2312 or newer, commSplit override**: set `EOF_COMMSPLIT=1` **before**
  `eofCompile` and keep it set when running. `eofCompile` then builds
  `libs/commSplit` into `FOAM_USER_LIBBIN`, the solvers link against it, and
  either order on the `mpirun` line works:
```
export EOF_COMMSPLIT=1
eofCompile
wmake $EOF_HOME/solvers/mhdInterFoamv2606
```
  Then copy the test case out of the repository, prepare it and run it:
```
cp -r $EOF_HOME/tests/levitation2Dv2606 ~/runs/
cd ~/runs/levitation2Dv2606
ElmerGrid 2 2 meshElmer -metis 2
setFields
decomposePar
mpirun -np 2 mhdInterFoamv2606 -parallel : -np 2 ElmerSolver_mpi case.sif
```
  The variable is not set by default; the library and the launch order must
  agree, so do not change it between building and running.

* **Tests**: `tests/travisTests_v2606.sh` builds the libraries and runs
  `interpolationTest` and `levitation2Dv2606`. It sources the OpenFOAM
  `etc/bashrc` given in `OPENFOAM_BASHRC` (default
  `/usr/lib/openfoam/openfoam2606/etc/bashrc`) and honours `EOF_COMMSPLIT`:
```
OPENFOAM_BASHRC=/path/to/OpenFOAM-v2606/etc/bashrc tests/travisTests_v2606.sh
```
  `tests/travisTests.sh` is the openfoam.org driver and skips the v2606 tests.
  The reference norm of `levitation2Dv2606` was produced with Elmer 9.0; newer
  Elmer releases may miss its tolerance slightly.
