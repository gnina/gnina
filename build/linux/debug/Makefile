BASE=/usr/local
NVCC=/usr/local/cuda/bin/nvcc 
LDFLAGS='-L/usr/lib -L. -L/usr/local/lib -L/usr/local/cuda/lib64'
C_PLATFORM= -pthread 
GPP=g++
C_OPTIONS= -O0  -g 
#C_OPTIONS= -O3 -DNDEBUG -g -flto
BOOST_LIB_VERSION=
NVOPTS=--expt-relaxed-constexpr -std=c++11 -gencode arch=compute_52,code=sm_52 -O0 -G  -Xcompiler -Wno-unreachable-code -gencode arch=compute_35,code=sm_35
OPENBABEL_INCLUDE=/usr/local/include/openbabel-2.0

include ../../makefile_common
