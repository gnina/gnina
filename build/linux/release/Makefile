BASE=/usr/local
NVCC=/usr/local/cuda/bin/nvcc 
LDFLAGS='-L/usr/lib -L. -L/usr/local/lib -L/usr/local/cuda/lib64'
C_PLATFORM= 
GPP=g++
C_OPTIONS= -O3 -DNDEBUG -g   
#C_OPTIONS= -O3 -DNDEBUG -g -flto
BOOST_LIB_VERSION=
NVOPTS=-gencode arch=compute_50,code=sm_50 -O3 -g -gencode arch=compute_35,code=sm_35
OPENBABEL_INCLUDE=/usr/local/include/openbabel-2.0
include ../../makefile_common
