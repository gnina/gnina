#pragma once 

#include <random>
#include <limits>
#include <algorithm>
#include <vector>
#include <functional>
#include "test/gnina/parsed_args.h"
#include "../loss.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

void test_cpu_l2();
void test_gpu_l2();
