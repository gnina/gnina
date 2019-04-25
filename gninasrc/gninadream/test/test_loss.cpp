#include "test_loss.h"
#define VECLEN 500

extern parsed_args p_args;

void test_cpu_l2() {
  p_args.log << "CPU L2 Test \n";
  p_args.log << "Using random seed: " << p_args.seed << '\n';
  //test 1: vector with itself
  std::mt19937 engine(p_args.seed);
  std::uniform_real_distribution<> vec_dist(-50, 50); 
  auto gen = std::bind(vec_dist, engine);
  std::vector<float> vec1(VECLEN);
  std::generate(vec1.begin(), vec1.end(), gen);
  float l2;
  cpu_l2(&vec1[0], &vec1[0], &l2, vec1.size());
  p_args.log << "true loss 0 calculated loss " << l2 << "\n";
  BOOST_CHECK_EQUAL(l2, 0);

  //test 2: distance between vec1 and vec2 = vec1 + 2
  std::vector<float> vec2(VECLEN);
  for (size_t i=0; i<VECLEN; ++i)
    vec2[i] = vec1[i] + 2;
  float dist = std::sqrt(4 * VECLEN);
  cpu_l2(&vec1[0], &vec2[0], &l2, vec1.size());
  p_args.log << "true loss " << dist << " calculated loss " << l2 << "\n";
  BOOST_REQUIRE_SMALL(l2 - dist, (float)0.001);
}

void test_gpu_l2() {
  p_args.log << "GPU L2 Test \n";
  p_args.log << "Using random seed: " << p_args.seed << '\n';
  //test 1: vector with itself
  std::mt19937 engine(p_args.seed);
  std::uniform_real_distribution<> vec_dist(-50, 50); 
  auto gen = std::bind(vec_dist, engine);
  std::vector<float> vec1(VECLEN);
  std::generate(vec1.begin(), vec1.end(), gen);
  float* gpu_v1;
  cudaMalloc(&gpu_v1, sizeof(float)*vec1.size());
  cudaMemcpy(gpu_v1, &vec1[0], sizeof(float)*vec1.size(),
      cudaMemcpyHostToDevice);
  float* gpu_l2;
  cudaMalloc(&gpu_l2, sizeof(float));
  do_gpu_l2(gpu_v1, gpu_v1, gpu_l2, vec1.size());
  float l2;
  cudaMemcpy(&l2, gpu_l2, sizeof(float), cudaMemcpyDeviceToHost);
  p_args.log << "true loss 0 calculated loss " << l2 << "\n";
  BOOST_CHECK_EQUAL(l2, 0);

  //test 2: distance between vec1 and vec2 = vec1 + 2
  std::vector<float> vec2(VECLEN);
  for (size_t i=0; i<VECLEN; ++i)
    vec2[i] = vec1[i] + 2;
  float* gpu_v2;
  cudaMalloc(&gpu_v2, sizeof(float)*vec2.size());
  cudaMemcpy(gpu_v2, &vec2[0], sizeof(float)*vec2.size(),
      cudaMemcpyHostToDevice);
  float dist = std::sqrt(4 * VECLEN);
  do_gpu_l2(gpu_v1, gpu_v2, gpu_l2, vec1.size());
  cudaMemcpy(&l2, gpu_l2, sizeof(float), cudaMemcpyDeviceToHost);
  p_args.log << "true loss " << dist << " calculated loss " << l2 << "\n";
  BOOST_REQUIRE_SMALL(l2 - dist, (float)0.001);
}
