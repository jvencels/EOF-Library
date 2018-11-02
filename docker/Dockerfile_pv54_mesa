FROM ubuntu:16.04

USER root

RUN apt-get -y update && \
    apt-get install -y   \
        pkg-config       \
        zlib1g-dev       \
        libexpat1-dev    \
        cmake            \
        mesa-utils       \
        python           \
        python-pip       \
        wget             \
        libopenmpi-dev   \
        openmpi-bin

RUN wget http://releases.llvm.org/7.0.0/llvm-7.0.0.src.tar.xz && \
    tar -xvf llvm-7.0.0.src.tar.xz && \
    cd llvm-7.0.0.src && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DLLVM_BUILD_LLVM_DYLIB=ON -DLLVM_ENABLE_RTTI=ON -DLLVM_INSTALL_UTILS=ON -DLLVM_TARGETS_TO_BUILD:STRING=X86 -DLLVM_LINK_LLVM_DYLIB=ON .. && \
    make -j4 install && \
    cd / && rm llvm-7.0.0.src.tar.xz && rm -r llvm-7.0.0.src

RUN wget https://mesa.freedesktop.org/archive/mesa-18.2.4.tar.gz && \
    tar -xvzf mesa-18.2.4.tar.gz && \
    cd mesa-18.2.4 && \
    ./configure --prefix=/usr/local/ --enable-opengl --disable-gles1 --disable-gles2 --disable-va --disable-xvmc --disable-vdpau --enable-shared-glapi --disable-texture-float --enable-llvm --enable-llvm-shared-libs --with-gallium-drivers=swrast --disable-dri --disable-glx --disable-egl --with-platforms= --disable-gbm --with-swr-archs=avx,avx2,skx,knl --disable-osmesa --enable-gallium-osmesa && \
    make -j4 install && \
    cd / && rm mesa-18.2.4.tar.gz && rm -r mesa-18.2.4

RUN wget "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.4&type=source&os=Sources&downloadFile=ParaView-v5.4.0.tar.gz" -O ParaView-v5.4.0.tar.gz && \
    tar -xvzf ParaView-v5.4.0.tar.gz && cd ParaView-v5.4.0 && mkdir build && cd build && \
    cmake -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_BUILD_QT_GUI=OFF -DVTK_USE_X=OFF -DOPENGL_INCLUDE_DIR=IGNORE -DOPENGL_xmesa_INCLUDE_DIR=IGNORE -DOPENGL_gl_LIBRARY=IGNORE -DOSMESA_INCLUDE_DIR=/usr/local/include -DOSMESA_LIBRARY=/usr/local/lib/libOSMesa.so -DVTK_OPENGL_HAS_OSMESA=ON -DVTK_USE_OFFSCREEN=OFF -DPARAVIEW_USE_MPI=ON .. && \
    make -j4 install && \
    cd / && rm ParaView-v5.4.0.tar.gz && rm -r ParaView-v5.4.0

