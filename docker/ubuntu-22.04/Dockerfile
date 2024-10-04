# Recommended build process

#get CUDA - getting an older version for compatibility, recommend changing to newer
FROM nvidia/cuda:12.0.0-devel-ubuntu22.04

# get packages 
WORKDIR /root
ARG DEBIAN_FRONTEND="noninteractive"
RUN apt-get update -y
RUN apt-get install build-essential -y wget -y libboost-all-dev -y libeigen3-dev -y libgoogle-glog-dev -y libprotobuf-dev -y protobuf-compiler -y libhdf5-dev -y libatlas-base-dev -y python3-dev -y librdkit-dev -y python3-numpy -y python3-pip -y python3-pytest -y swig -y
RUN apt update ; apt upgrade -y ; apt-get -y install build-essential git wget libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest libjsoncpp-dev
RUN apt-get update -y && apt-get install git-all -y curl && apt-get update -y
RUN pip3 install cmake scikit-image pyquaternion google-api-python-client six
RUN pip3 install torch torchvision torchaudio

RUN git clone https://github.com/openbabel/openbabel.git && cd openbabel &&  git checkout openbabel-3-1-1 && mkdir build && cd build && cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON .. && make -j4 && make install

ADD "https://github.com/gnina/gnina/commits?per_page=1" latest_commit

RUN git clone https://github.com/gnina/gnina.git; \
    cd gnina; mkdir build; cd build; \
    cmake .. ;\
    make -j4 ; make install 
