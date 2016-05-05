BASE=/usr/local/
C_PLATFORM= -pthread  -DSMINA_GPU
GPP=g++
C_OPTIONS= -O0  -g   
BOOST_LIB_VERSION=
NVCC=nvcc -std=c++11 -Xcudafe "--diag_suppress=boolean_controlling_expr_is_constant"
NVOPTS=-gencode arch=compute_52,code=sm_52 -O0 -g
OPENBABEL_INCLUDE=/usr/include/openbabel-2.0
GPUOBJ = gpucode.o bfgs.o quasi_newton.o

include ../../makefile_common

.DEFAULT_GOAL := gpu

%.o : ../../../src/lib/%.cu 
	$(NVCC) $(NVOPTS) -DSMINA_GPU -x cu -dc $< -o $@

gpu: $(MAINOBJ) $(LIBOBJ) $(GPUOBJ)
	$(NVCC) $(NVOPTS) $(LDFLAGS) -o smina.$@ $^ $(LIBS)  -lcudart
