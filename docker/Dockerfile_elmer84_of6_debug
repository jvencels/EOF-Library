ARG BASE_IMAGE=hfdresearch/swak4foamandpyfoam:latest-v6.0

# Alternatively, you can use Docker flag
# --build-arg BASE_IMAGE=openfoam/openfoam6-paraview54:latest

FROM $BASE_IMAGE

USER root

RUN apt-get -y update && \
  apt-get install -y \
    gcc \
    cmake \
    g++ \
    gfortran \
    libopenblas-dev \
    git \
    valgrind

# OpenMPI
RUN cd /opt && \
  wget https://download.open-mpi.org/release/open-mpi/v1.10/openmpi-1.10.2.tar.gz && \
  tar xvzf openmpi-1.10.2.tar.gz && \
  rm openmpi-1.10.2.tar.gz && \
  cd openmpi-1.10.2 && \
  ./configure --prefix=$(pwd)/install && \
  make -j install && \
  echo "export MPI_INST=/opt/openmpi-1.10.2/install" >> /etc/bash.bashrc && \
  echo "export PATH=\$MPI_INST/bin:\$MPI_INST/include:\$PATH" >> /etc/bash.bashrc && \
  echo "export LD_LIBRARY_PATH=\$MPI_INST/lib:\$LD_LIBRARY_PATH" >> /etc/bash.bashrc

# Elmer
RUN cd /opt && \
  git clone https://github.com/ElmerCSC/elmerfem.git && \
  chown -R openfoam:openfoam elmerfem && \
  cd elmerfem && \
  git checkout release && \
  git log --pretty=oneline | head -n 10 | tee -a /opt/elmerBuild.log && \
  mkdir build && \
  cd build && \
  cmake -DWITH_MPI=TRUE -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS_DEBUG="-cpp -g -O0 -fbacktrace -fcheck=all -Wall -Wtabs -Wuninitialized -Wextra -Wconversion -ffpe-trap=overflow -ffree-line-length-none" ../ | tee -a /opt/elmerBuild.log && \
  make -j install | tee -a /opt/elmerBuild.log
