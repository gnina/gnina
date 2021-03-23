# Recommended build process

#get CUDA and cudnn (cudnn8 currently not supported)
FROM nvidia/cuda:10.2-cudnn7-devel-ubuntu18.04

# get packages 
WORKDIR /root
RUN apt update ; apt upgrade -y ; apt-get -y install build-essential git wget libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest


# set up acceptable cmake version
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.6/cmake-3.18.6-Linux-x86_64.tar.gz; \
    tar xzf cmake-3.18.6-Linux-x86_64.tar.gz; \
    cp cmake-3.18.6-Linux-x86_64/bin/* /usr/local/bin; \
    cp -r cmake-3.18.6-Linux-x86_64/share/cmake-3.18 /usr/local/share/

# build and install openbabel 3
RUN git clone https://github.com/openbabel/openbabel.git; \
    cd openbabel; git checkout openbabel-3-1-1; mkdir build; cd build;  \
    cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF ..;  \
    make -j4; make install


RUN git clone https://github.com/gnina/gnina.git; \
    cd gnina; mkdir build; cd build; \
    cmake .. ;\
    make -j4 ; make install 
