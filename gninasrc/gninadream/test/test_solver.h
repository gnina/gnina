#pragma once 
#include "caffe/sgd_solvers.hpp"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

void test_iopt_update(caffe::Solver<float>* solver);

void test_iopt_improvement(caffe::Solver<float>* solver);

void test_iopt_exclude_rec(caffe::Solver<float>* solver);
