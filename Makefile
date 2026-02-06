export ELMER_HOME        := /usr/local
export ELMER_SOLVER_HOME := /usr/local
export EOF_HOME          := /home/openfoam/EOF-Library
export EOF_SRC           := $(EOF_HOME)/libs
export PATH              := /usr/local/bin:$(PATH)
export LD_LIBRARY_PATH   := /usr/local/lib:$(LD_LIBRARY_PATH)
export OPENFOAM_HOME	 := /opt/openfoam6
export ROOT_DIR			 := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST)))) # https://stackoverflow.com/a/23324703
export IMAGE_NAME		 := eof_local


# This is so that the environment variables persist between commands
SHELL := /bin/bash
.ONESHELL:

# -- Docker Container

environment:
	. $(OPENFOAM_HOME)/etc/bashrc
	. $(EOF_HOME)/etc/bashrc
	cd $(EOF_HOME)

eof: environment
	. $(OPENFOAM_HOME)/etc/bashrc && wclean $(EOF_SRC)/coupleElmer
	. $(OPENFOAM_HOME)/etc/bashrc && wmake $(EOF_SRC)/coupleElmer
	elmerf90 -o $(EOF_SRC)/Elmer2OpenFOAM.so -J $(nproc) $(EOF_SRC) $(EOF_SRC)/Elmer2OpenFOAM.F90
	elmerf90 -o $(EOF_SRC)/OpenFOAM2Elmer.so -J $(nproc) $(EOF_SRC) $(EOF_SRC)/OpenFOAM2Elmer.F90

solver: environment
	. $(OPENFOAM_HOME)/etc/bashrc && wclean solvers/mdhLinearHall
	. $(OPENFOAM_HOME)/etc/bashrc && wmake solvers/mdhLinearHall
	rm -rf solvers/mdhLinearHall/processor*
	echo "FOAM_USER_APPBIN=$FOAM_USER_APPBIN" >> log.txt
	echo "FOAM_APPBIN=$FOAM_APPBIN" >> log.txt

simulation: environment 
	cd $(EOF_HOME)/tests/linearHall
	rm -rf processor* 
	. $(EOF_HOME)/etc/bashrc && . $(OPENFOAM_HOME)/etc/bashrc && blockMesh && potentialFoam && decomposePar
	ElmerGrid 14 2 channel.msh -out meshElmer -autoclean -merge 1e-8 -removeunused
	ElmerGrid 2 2 meshElmer -metis 2
	echo \"case.sif\" > ELMERSOLVER_STARTINFO
	cd $(EOF_HOME) 

elmer: environment
	sudo cp libs/solvers/MHDSolve/MHDSolve.F90 /opt/elmerfem/fem/src/modules/MHDSolve.F90 
	cd /opt/elmerfem/build
	sudo cmake .. 
	sudo make MHDSolve 
	sudo make install MHDSolve 
	cd $(EOF_HOME)

# -- Host System

build_environment:
	cd $(ROOT_DIR)

setup: build_environment
	mkdir -p ./runs

build: setup
	docker build --progress=plain --network host --platform linux/amd64 -f docker/Dockerfile.build_simulation -t $(IMAGE_NAME):latest .

clean: build_environment
	docker ps -a --filter "ancestor=$(IMAGE_NAME)" -q | xargs -r docker rm -f
	docker images $(IMAGE_NAME) -q | xargs -r docker rmi -f
	rm -rf runs
