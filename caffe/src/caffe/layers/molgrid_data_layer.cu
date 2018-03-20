#include <stdint.h>
#include <vector>


#include "caffe/layers/molgrid_data_layer.hpp"


//gridding is implemented in gridmaker
#include "gninasrc/lib/gridmaker.cu"

namespace caffe {

template <typename Dtype>
void MolGridDataLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
		const vector<Blob<Dtype>*>& top)
{
	forward(bottom, top, true);
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
		const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
	backward(top, bottom, true);
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::setAtomGradientsGPU(GridMaker& gmaker, Dtype *diff)  {

  //launch a kernel for each batch
  for (int item_id = 0; item_id < batch_size; ++item_id) {
    int offset = item_id*example_size;
    //malloc and copy batch data
    float4* atoms;
    short* whichGrid;
    float3* gradient;

    int natoms = transform.mol.atoms.size();
    mol_transform& transform = batch_transform[item_id];
    cudaMalloc(&atoms, sizeof(float4)*natoms);
    cudaMemcpy(atoms, &transform.mol.atoms[0],
            sizeof(float4)*transform.mol.atoms.size(), cudaMemcpyHostToDevice);
    cudaMalloc(&whichGrid, sizeof(short)*transform.mol.whichGrid.size());
    cudaMemcpy(whichGrid, &transform.mol.whichGrid[0],
            sizeof(short)*transform.mol.whichGrid.size(),
            cudaMemcpyHostToDevice);
    cudaMalloc(&gradient, sizeof(float3)*transform.mol.gradient.size());
    cudaMemset(gradient, 0, sizeof(float3)*transform.mol.gradient.size());
    gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);

    //diff is batch x channel x X x Y x Z
    setAtomGradientsGPU<<<1, natoms>>>(gmaker, atoms, whichGrid, gradient, 
            transform.mol.center, transform.mol.Q, transform.center, diff,
            offset);
  }
}

INSTANTIATE_LAYER_GPU_FORWARD(MolGridDataLayer);
INSTANTIATE_LAYER_GPU_BACKWARD(MolGridDataLayer);

}  // namespace caffe
