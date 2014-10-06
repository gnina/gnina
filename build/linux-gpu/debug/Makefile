BASE=/usr/local/
C_PLATFORM= -pthread  -DSMINA_GPU
GPP=g++
C_OPTIONS= -O0  -g   
BOOST_LIB_VERSION=
NVCC=nvcc
NVOPTS=-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -O0 -g
OPENBABEL_INCLUDE=/usr/include/openbabel-2.0
GPUOBJ = gpucode.o

include ../../makefile_common

.DEFAULT_GOAL := gpu

%.o : ../../../src/lib/%.cu 
	$(NVCC) $(NVOPTS) -DSMINA_GPU -o $@ -c $< 
	
gpu: $(MAINOBJ) $(LIBOBJ) $(GPUOBJ)
	$(CC) $(LDFLAGS) -o smina.$@ $^ $(LIBS)  -lcudart
