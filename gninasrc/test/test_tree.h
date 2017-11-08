#pragma once

void make_tree(model* m);

__global__ 
void increment_kernel(conf_gpu c, const change_gpu g, fl factor, gpu_data* gdata);

__global__ 
void set_conf_kernel(gpu_data* gdata, const conf_gpu c);

void test_set_conf();
