BASE=/usr/local
NVCC=/usr/local/cuda/bin/nvcc 
LDFLAGS='-L/usr/lib -L. -L/usr/local/lib -L/usr/local/cuda/lib64'
C_PLATFORM= -pthread 
GPP=g++
C_OPTIONS= -O3 -DNDEBUG -g   
#C_OPTIONS= -O3 -DNDEBUG -g -flto
BOOST_LIB_VERSION=
NVOPTS=-std=c++11 -gencode arch=compute_52,code=sm_52 -O3
OPENBABEL_INCLUDE=/usr/local/include/openbabel-2.0

include ../../makefile_common

.DEFAULT_GOAL := gpu

%.o : ../../../src/lib/%.cu 
		$(NVCC) $(NVOPTS) -DSMINA_GPU -o $@ -dc $<

gpu: $(MAINOBJ) $(LIBOBJ) $(GPUOBJ)
		$(NVCC) -gencode arch=compute_52,code=sm_52 $(LDFLAGS) -o smina.$@ $^ $(LIBS)  -lcudart

