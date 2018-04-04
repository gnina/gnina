#include "test_utils.h"
#include "test_cnn.h"
#include "atom_constants.h"
#include "gridmaker.h"
#include "cnn_scorer.h"
#include <cuda_runtime.h>
#include <boost/math/quaternion.hpp>
#include "quaternion.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "caffe/util/rng.hpp"

extern parsed_args p_args;

typedef boost::math::quaternion<float> quaternion;
void test_set_atom_gradients() {
    // randomly generate gridpoint gradients, accumulate for atoms, and then compare
    p_args.log << "CNN Set Atom Gradients Test \n";
    p_args.log << "Using random seed: " << p_args.seed << '\n';
    p_args.log << "Iteration " << p_args.iter_count << '\n';
    std::mt19937 engine(p_args.seed);
    std::normal_distribution<> diff_dist(0, 5);
    auto gen = std::bind(diff_dist, engine);

    // TODO? not using the actual CNN atom types, but I don't think this should
    // matter for correctness. however, I did template the relevant test_utils
    // functions to facilitate changing this. the following also (relatedly) 
    // does not distinguish between ligand and receptor atoms.
    std::vector<atom_params> mol_atoms;
    std::vector<smt> mol_types;
    make_mol(mol_atoms, mol_types, engine, 0, 1, 5000, 11.5, 11.5, 11.5);
    cnn_options cnnopts;
    cnnopts.cnn_scoring = true;
    model m;
    CNNScorer cnn_scorer(cnnopts, m);
    caffe::MolGridDataLayer<CNNScorer::Dtype>* mgrid = cnn_scorer.mgrid;
    mgrid->batch_transform.resize(1);
    mgrid->batch_transform[0] = caffe::MolGridDataLayer<CNNScorer::Dtype>::mol_transform();
    caffe::MolGridDataLayer<CNNScorer::Dtype>::mol_transform& transform = mgrid->batch_transform[0];
    vec center(0,0,0);
    for (size_t i=0; i<mol_atoms.size(); ++i) {
        atom_params& ainfo = mol_atoms[i];
        transform.mol.atoms.push_back(make_float4(ainfo.coords.x,
                    ainfo.coords.y, ainfo.coords.z, xs_radius(mol_types[i])));
        transform.mol.whichGrid.push_back(mol_types[i]);
        transform.mol.gradient.push_back(make_float3(0,0,0));
        center += vec(ainfo.coords.x, ainfo.coords.y, ainfo.coords.z);
    }
    center /= mol_atoms.size();
    transform.center = center;
    transform.Q = quaternion(1,0,0,0);
    transform.set_random_quaternion(caffe::caffe_rng());
    //initialize gmaker
    GridMaker gmaker;
    bool spherize = false;
    double dim = round(mgrid->dimension/mgrid->resolution)+1;
    gmaker.initialize(mgrid->resolution, mgrid->dimension, mgrid->radiusmultiple, 
            mgrid->binary, spherize);
    gmaker.setCenter(center[0], center[1], center[2]);

    //randomly intialize gridpoint gradients
    std::vector<float> diff((smt::NumTypes)*dim*dim*dim);
    generate(begin(diff), end(diff), gen);

    //set up and calculate CPU atom gradients
    caffe::MolGridDataLayer<CNNScorer::Dtype>::Grids grids(&diff[0], 
            boost::extents[smt::NumTypes][dim][dim][dim]);
    caffe::MolGridDataLayer<CNNScorer::Dtype>::mol_transform cpu_transform =
        mgrid->batch_transform[0];
    gmaker.setAtomGradientsCPU(cpu_transform.mol.atoms, cpu_transform.mol.whichGrid,
            cpu_transform.Q, grids, cpu_transform.mol.gradient);

    //calculate GPU atom gradients
    float* gpu_diff;
    cudaMalloc(&gpu_diff, diff.size() * sizeof(float));
    cudaMemcpy(gpu_diff, &diff[0], diff.size() * sizeof(float),
            cudaMemcpyHostToDevice);
    mgrid->setAtomGradientsGPU(gmaker, gpu_diff, 1);

    //compare results
    for (size_t i=0; i<mol_atoms.size(); ++i) {
        for (size_t j=0; j<3; ++j) {
            p_args.log << "CPU " << cpu_transform.mol.gradient[i][j] << 
                " GPU " << transform.mol.gradient[i][j] << "\n";
            BOOST_REQUIRE_SMALL(cpu_transform.mol.gradient[i][j]- 
                    transform.mol.gradient[i][j], (float)0.01);
        }
    }
    cudaFree(gpu_diff);
}
