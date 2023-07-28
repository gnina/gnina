# Recommended build process

#get CUDA 
FROM nvidia/cuda:11.7.0-devel-ubuntu22.04

# get packages 
WORKDIR /root
RUN apt update ; apt upgrade -y ; apt-get -y install build-essential git cmake wget libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest openbabel libopenbabel-dev python3-openbabel

ADD "https://github.com/gnina/gnina/commits?per_page=1" latest_commit

RUN git clone https://github.com/gnina/gnina.git; \
    cd gnina; mkdir build; cd build; \
    cmake -DCUDA_ARCH_NAME=All .. ;\
    make -j4 ; make install 
