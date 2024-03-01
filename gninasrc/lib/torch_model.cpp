/*
 * torch_model.cpp
 *
 *  Created on: Feb 14, 2024
 *      Author: dkoes
 */

#include "torch_model.h"
#include "common.h"
#include <json/json.h>
#include <string>

using namespace std;
using namespace libmolgrid;

static const string default_recmap(R"(AliphaticCarbonXSHydrophobe 
AliphaticCarbonXSNonHydrophobe 
AromaticCarbonXSHydrophobe 
AromaticCarbonXSNonHydrophobe
Bromine Iodine Chlorine Fluorine
Nitrogen NitrogenXSAcceptor 
NitrogenXSDonor NitrogenXSDonorAcceptor
Oxygen OxygenXSAcceptor 
OxygenXSDonorAcceptor OxygenXSDonor
Sulfur SulfurAcceptor
Phosphorus 
Calcium
Zinc
GenericMetal Boron Manganese Magnesium Iron
)");

static const string default_ligmap(R"(AliphaticCarbonXSHydrophobe 
AliphaticCarbonXSNonHydrophobe 
AromaticCarbonXSHydrophobe 
AromaticCarbonXSNonHydrophobe
Bromine Iodine
Chlorine
Fluorine
Nitrogen NitrogenXSAcceptor 
NitrogenXSDonor NitrogenXSDonorAcceptor
Oxygen OxygenXSAcceptor 
OxygenXSDonorAcceptor OxygenXSDonor
Sulfur SulfurAcceptor
Phosphorus
GenericMetal Boron Manganese Magnesium Zinc Calcium Iron
)");

/** Read in a torch script model from the provided stream */
template <bool isCUDA> TorchModel<isCUDA>::TorchModel(std::istream &in, const string &name, tee *log) {
  c10::Device device = isCUDA ? torch::kCUDA : torch::kCPU;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    torch::jit::ExtraFilesMap extras;
    extras["metadata"] = ""; // only loads if key already present
    module = torch::jit::load(in, device, extras);
    module.to(device);

    string data = extras["metadata"];
    string recmap, ligmap;
    double resolution = 0.5;
    double dimension = 23.5;
    if (data.length() == 0) {
      if (log)
        *log << "WARNING: Torch model missing metadata.  Default grid parameters will be used.\n";
      recmap = default_recmap;
      ligmap = default_ligmap;
    } else {
      Json::Reader reader;
      Json::Value root;
      reader.parse(data, root);

      if (root.isMember("resolution")) {
        resolution = root["resolution"].asDouble();
      } else {
        if (log)
          *log << "WARNING: Resolution not specified in model file.  Using default.\n";
      }
      if (root.isMember("dimension")) {
        dimension = root["dimension"].asDouble();
      } else {
        if (log)
          *log << "WARNING: Dimension not specified in model file.  Using default.\n";
      }
      if (root.isMember("ligmap")) {
        ligmap = root["ligmap"].asString();
      } else {
        if (log)
          *log << "WARNING: ligmap not specified in model file.  Using default.\n";
      }
      if (root.isMember("recmap")) {
        recmap = root["recmap"].asString();
      } else {
        if (log)
          *log << "WARNING: recmap not specified in model file.  Using default.\n";
      }
      if (root.isMember("skip_loss")) {
        skip_loss = root["skip_loss"].asBool();
      }
    }

    gmaker.initialize(resolution, dimension);

    // setup typers
    stringstream ligstream(ligmap), recstream(recmap);
    lig_typer = make_shared<FileMappedGninaTyper>(ligstream);
    rec_typer = make_shared<FileMappedGninaTyper>(recstream);

  } catch (const c10::Error &e) {
    throw usage_error("Could not read torch model " + name);
  }
}

static CoordinateSet make_coordset(const vector<float3> &coords, const vector<smt> &smtypes,
                                   shared_ptr<AtomTyper> typer) {
  if (coords.size() != smtypes.size())
    throw internal_error("Shape mismatch", __LINE__);

  vector<float> types;
  types.reserve(coords.size());
  vector<float> radii;
  radii.reserve(coords.size());
  for (unsigned i = 0, n = smtypes.size(); i < n; i++) {
    smt origt = smtypes[i];
    auto t_r = typer->get_int_type(origt);
    int t = t_r.first;
    types.push_back(t);
    radii.push_back(t_r.second);

    if (t < 0 && origt > 1) { // don't warn about hydrogens
      std::cerr << "Unsupported ligand atom type " << GninaIndexTyper::gnina_type_name(origt) << "\n";
    }
  }

  return CoordinateSet(coords, types, radii, typer->num_types());
}

//wrapper to get appropriate grid from an MGrid for template value of isCUDA
template <bool isCUDA>
static Grid<float, 2, isCUDA> get2DGrid(MGrid2f& mgrid);
template<>
Grid<float, 2, true> get2DGrid(MGrid2f& mgrid) { return mgrid.gpu();}
template<>
Grid<float, 2, false> get2DGrid(MGrid2f& mgrid) { return mgrid.cpu();}

template <bool isCUDA>
std::vector<float> TorchModel<isCUDA>::forward(const std::vector<float3> &rec_coords, const std::vector<smt> &rec_types,
                                               const std::vector<float3> &lig_coords, const std::vector<smt> &lig_types,
                                               const vec &center, bool rotate, bool compute_gradient) {

  torch::AutoGradMode enable_grad(compute_gradient);
  // make coordinate sets
  CoordinateSet rec = make_coordset(rec_coords, rec_types, rec_typer);
  CoordinateSet lig = make_coordset(lig_coords, lig_types, lig_typer);

  // set center from ligand if not specified
  float3 gcenter = {center.x(), center.y(), center.z()};
  if (!isfinite(center.x())) {
    gcenter = lig.center();
  }

  CoordinateSet combined(rec, lig);

  Transform transform(gcenter, 0, rotate);
  if (rotate) {
    transform.forward(combined, combined);
  }
  // create grids
  long ntypes = combined.num_types();
  long gd = gmaker.get_first_dim();

  auto options = torch::TensorOptions().dtype(torch::kFloat32).device(isCUDA ? torch::kCUDA : torch::kCPU).requires_grad(compute_gradient);
  torch::Tensor gtensor = torch::zeros({1, ntypes, gd, gd, gd}, options);
  Grid<float, 4, isCUDA> out(gtensor.data_ptr<float>(), ntypes, gd, gd, gd);
  gmaker.forward(gcenter, combined, out);

  // evaluate model
  vector<torch::jit::IValue> inputs{gtensor};
  auto result = module.forward(inputs).toTuple()->elements();

  // get results
  auto pose_logit = result[0].toTensor();
  auto pose = torch::softmax(pose_logit, 1).index({0, 1}).item<float>();
  auto affinity = result[1].toTensor()[0].item<float>();

  auto loptions = torch::TensorOptions().dtype(torch::kLong).device(isCUDA ? torch::kCUDA : torch::kCPU);
  torch::Tensor labels = torch::ones({1}, loptions);

  auto loss = skip_loss ? pose_logit.index({0,1}) : torch::cross_entropy_loss(pose_logit, labels);

  if (compute_gradient) {
    loss.backward(); 
    auto grad = gtensor.grad();
    Grid<float, 4, isCUDA> gridgrad(grad.data_ptr<float>(), ntypes, gd, gd, gd);
    MGrid2f atomic_gradients(combined.size(),3);
    auto coord_grad = get2DGrid<isCUDA>(atomic_gradients);
    gmaker.backward(gcenter,combined, gridgrad, coord_grad);
    if(rotate) {
      transform.backward(coord_grad,coord_grad,false);
    }
    unsigned nr = rec.size();
    unsigned nl = lig.size();
    VINA_CHECK(nr+nl == combined.size());
    gradient_rec.resize(nr);

    auto cg_cpu = atomic_gradients.cpu();
    for(unsigned i = 0; i < nr; i++) {
      gradient_rec[i] = gfloat3(cg_cpu[i][0],cg_cpu[i][1],cg_cpu[i][2]);
    }
    gradient_lig.resize(nl);
    for(unsigned i = 0; i < nl; i++) {
      gradient_lig[i] = gfloat3(cg_cpu[i+nr][0],cg_cpu[i+nr][1],cg_cpu[i+nr][2]);
    }

  }
  vector<float> scores{pose, affinity, loss.item<float>()};
  return scores;
}

template <bool isCUDA> void TorchModel<isCUDA>::getLigandGradient(std::vector<gfloat3> &grad) {
  grad = gradient_lig;
}

template <bool isCUDA> void TorchModel<isCUDA>::getReceptorGradient(std::vector<gfloat3> &grad) {
  grad = gradient_rec;
}

// explicit instaniationsº
template class TorchModel<true>;
template class TorchModel<false>;
