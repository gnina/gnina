#ifndef CAFFE_MATMUL_LAYER_HPP_
#define CAFFE_MATMUL_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

namespace caffe {

/**
 * @brief Matrix multiplication of two bottom blobs.
 *
 * Computes a batch of matrix-matrix multiplications between two bottom
 * blobs after implicitly flattening each bottom blob after axis 1 and
 * optionally transposing them. The computation is as follows:
 *
 * C[...] = 0.0
 * for b in range(B):
 *     for m in range(M):
 *         for n in range(N):
 *             for k in range(K):
 *                 C[b,m,n] += A[b,m,k] * B[b,k,n]
 *
 * Axis 0, the batch size, must be the same for each bottom.
 * Axis 1 is the number of rows of the matrices.
 * Axis 2 and beyond are the number of columns of the matrices.
 * 
 * After transposing, if applicable, the number of columns in the
 * first matrix must equal the number of rows in the second.
 * Then the top blob shape will have three axes corresponding to
 * the batch size, the number of rows in the first (transposed) matrix,
 * and the number of columns in the second (transposed) matrix.
 */
template <typename Dtype>
class MatMulLayer : public Layer<Dtype> {
 public:
  explicit MatMulLayer(const LayerParameter& param)
      : Layer<Dtype>(param) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "MatMul"; }
  virtual inline int ExactNumBottomBlobs() const { return 2; }
  virtual inline int ExactNumTopBlobs() const { return 1; }

 protected:
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_relevance(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  int batch_size; // num of matrix pairs to multiply
  int M; // num of rows of first matrix
  int K; // num of cols of first matrix/rows of second matrix
  int N; // num cols of second matrix

  bool transpose_A; // whether to transpose first input
  bool transpose_B; // whether to transpose second input
};

}  // namespace caffe

#endif  // CAFFE_MATMUL_LAYER_HPP_
