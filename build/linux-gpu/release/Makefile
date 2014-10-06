BASE=/usr/local/
C_PLATFORM= -pthread  -DSMINA_GPU
GPP=g++
C_OPTIONS= -O3 -DNDEBUG -g   
#C_OPTIONS= -O3 -DNDEBUG -g -flto
BOOST_LIB_VERSION=
NVCC=nvcc
NVOPTS=-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -O3
OPENBABEL_INCLUDE=/usr/include/openbabel-2.0
GPUOBJ = gpucode.o 

include ../../makefile_common

.DEFAULT_GOAL := gpu

%.o : ../../../src/lib/%.cu ../../../src/lib/*gpu.h
	$(NVCC) $(NVOPTS) -DSMINA_GPU -o $@ -c $< 
	
gpu: $(MAINOBJ) $(LIBOBJ) $(GPUOBJ)
	$(CC) $(LDFLAGS) -o smina.$@ $^ $(LIBS)  -lcudart
