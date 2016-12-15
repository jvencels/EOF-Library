# OpenFOAM_Elmer
Modiefied OpenFOAM library and test cases for MPI-coupled OpenFOAM + Elmer solvers. Used for solving coupled electromagnetic induction and MHD problems.

## Requirements ##

Tested on:
* Ubuntu 14.04
* Open MPI v2.0.1 (https://www.open-mpi.org/software/ompi/v2.0/)

## How to ##
* Download Elmer (https://github.com/ElmerCSC/elmerfem.git)
* Change git branch to "ElmerCommWorld"
* Obtain `Elmer2OpenFOAM.F90` solver (currently not publicly published)

```
git checkout ElmerCommWorld
 ```

* Compile Elmer with `-DWITH_MPI=TRUE`
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
cd ../mhdPisoFoam
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
cp -r ../tests/mhdPisoFoamTest .
```

* Prepare test

```
cd mhdPisoFoamTest
ElmerGrid 8 2 mesh_Elmer
ideasUnvToFoam mesh_OpenFOAM.unv
decomposePar
```

* Run OpenFOAM on 2 processes and Elmer on 1 process

```
mpirun -np 2 mhdPisoFoam -parallel : -np 1 ElmerSolver_mpi
```

* Postprocessing

```
reconstructPar
paraFoam
```

* Tick "JxB" and "JH" fields and load OpenFOAM results
* To open Elmer results, open `mesh_Elmer/case_0001.vtu`
