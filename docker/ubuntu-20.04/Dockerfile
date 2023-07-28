# Recommended build process

#get CUDA 
FROM nvidia/cuda:12.2.0-devel-ubuntu20.04

# get packages 

WORKDIR /root
RUN apt update
# https://stackoverflow.com/questions/44331836/apt-get-install-tzdata-noninteractive
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata

RUN apt update ; apt upgrade -y ; apt-get -y install build-essential git cmake wget libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest openbabel python3-openbabel

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
    
ADD "https://github.com/gnina/gnina/commits?per_page=1" latest_commit

RUN git clone https://github.com/gnina/gnina.git; \
    cd gnina; mkdir build; cd build; \
    cmake -DCUDA_ARCH_NAME=All .. ;\
    make -j4 ; make install 
