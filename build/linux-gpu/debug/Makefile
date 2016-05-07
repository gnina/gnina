BASE=/usr/local/
C_PLATFORM= -pthread  -DSMINA_GPU
GPP=g++
C_OPTIONS= -O0  -g   
BOOST_LIB_VERSION=
NVCC=nvcc
NVOPTS=-std=c++11 -gencode arch=compute_52,code=sm_52 -O0 -g
OPENBABEL_INCLUDE=/usr/include/openbabel-2.0
GPUOBJ = gpucode.o

include ../../makefile_common

.DEFAULT_GOAL := gpu

%.o : ../../../src/lib/%.cu 
	$(NVCC) $(NVOPTS) -DSMINA_GPU -o $@ -dc $<

gpu: $(MAINOBJ) $(LIBOBJ) $(GPUOBJ)
	$(NVCC) -gencode arch=compute_52,code=sm_52 $(LDFLAGS) -o smina.$@ $^ $(LIBS)  -lcudart
