#include "dl_scorer.h"
#include "caffe/caffe.hpp"

// set the default device to device and return cuda error code if there's a problem
int initializeCUDA(int device) {
  cudaError_t error;
  cudaDeviceProp deviceProp;

  error = cudaSetDevice(device);
  if (error != cudaSuccess) {
    // be silent if GPU not present
    return error;
  }

  error = cudaGetDevice(&device);

  if (error != cudaSuccess) {
    std::cerr << "cudaGetDevice returned error code " << error << "\n";
    return error;
  }

  error = cudaGetDeviceProperties(&deviceProp, device);

  if (deviceProp.computeMode == cudaComputeModeProhibited) {
    std::cerr << "Error: device is running in <Compute Mode Prohibited>, no threads can use ::cudaSetDevice().\n";
    return -1;
  }

  if (error != cudaSuccess) {
    return error;
  }

  caffe::Caffe::SetDevice(device);
  caffe::Caffe::set_mode(caffe::Caffe::GPU);
  return 0;
}

// Extract ligand atoms and coordinates
void DLScorer::setLigand(const model &m) {
  const atomv &atoms = m.get_movable_atoms();
  const vecv &coords = m.coordinates();
  num_atoms = atoms.size();

  if (m.ligands.size() == 0 && m.coordinates().size() == 0) {
    return;
  } else if (m.ligands.size() == 0) { // check for covalent
    ligand_smtypes.resize(0);
    ligand_coords.resize(0);
    ligand_map.resize(0);

    sz num_flex_atoms = m.m_num_movable_atoms;
    if (m.ligands.size() > 0) {
      num_flex_atoms = m.ligands[0].node.begin;
    }

    for (sz i = 0; i < num_flex_atoms; i++) {
      if (atoms[i].iscov) {
        const vec &c = coords[i];
        ligand_smtypes.push_back(atoms[i].sm);
        ligand_coords.push_back(float3({c[0], c[1], c[2]}));
        ligand_map.push_back(i);
      }
    }

    // Append inflex
    for (sz i = m.m_num_movable_atoms, n = coords.size(); i < n; i++) {
      if (atoms[i].iscov) {
        const vec &c = coords[i];
        ligand_smtypes.push_back(atoms[i].sm);
        ligand_coords.push_back(float3({c[0], c[1], c[2]}));
        ligand_map.push_back(i);
      }
    }
  } else {
    // Ligand atoms start at ligand root node start idx
    sz n = m.m_num_movable_atoms - m.ligands[0].node.begin;

    // Get ligand types and radii
    ligand_smtypes.resize(n);
    ligand_coords.resize(n);
    ligand_map.resize(n);

    sz offset = m.ligands[0].node.begin;
    for (sz i = 0; i < n; ++i) {
      ligand_smtypes[i] = atoms[i + offset].sm;
      const vec &coord = coords[i + offset];
      ligand_map[i] = i + offset;
      ligand_coords[i] = gfloat3(coord[0], coord[1], coord[2]);
    }
  }
}

// Extracts receptor atoms and coordinates
// Flex and inflex coordinates are taken from the model's movable atoms
// Flex coordinates are stored at the beginning, then inflex, then fixed
void DLScorer::setReceptor(const model &m) {
  num_atoms = m.get_movable_atoms().size();
  auto cbegin = m.get_movable_atoms().cbegin();
  auto flexend = cbegin + m.m_num_movable_atoms;
  if (m.ligands.size() > 0) {
    flexend = cbegin + m.ligands[0].node.begin;
  }
  // Number of receptor movable atoms
  sz num_flex_atoms = std::distance(cbegin, flexend);

  // Number of inflex atoms
  sz n_inflex = std::distance(m.get_movable_atoms().cbegin() + m.m_num_movable_atoms, m.get_movable_atoms().cend());

  // Number of fixed receptor atoms
  sz n_rigid = m.get_fixed_atoms().size();

  // Total receptor size
  sz n = num_flex_atoms + n_inflex + n_rigid;
  sz ncov = 0;
  if (receptor_smtypes.size() == 0) { // Do once at setup

    receptor_smtypes.reserve(n);
    // Insert flexible residues movable atoms
    for (auto it = cbegin; it != flexend; ++it) {
      smt origt = it->sm; // Original smina type
      if (!it->iscov) {
        receptor_smtypes.push_back(origt);
      } else {
        ncov += 1;
      }
    }

    CHECK_EQ(receptor_smtypes.size(), num_flex_atoms - ncov);

    // Insert inflex atoms
    cbegin = m.get_movable_atoms().cbegin() + m.m_num_movable_atoms;
    auto cend = m.get_movable_atoms().cend();
    for (auto it = cbegin; it != cend; ++it) {
      smt origt = it->sm; // Original smina type
      if (!it->iscov) {
        receptor_smtypes.push_back(origt);
      } else {
        ncov += 1;
      }
    }

    CHECK_EQ(receptor_smtypes.size(), num_flex_atoms + n_inflex - ncov);

    // Insert fixed receptor atoms
    cbegin = m.get_fixed_atoms().cbegin();
    cend = m.get_fixed_atoms().cend();
    for (auto it = cbegin; it != cend; ++it) {
      smt origt = it->sm; // Original smina type
      receptor_smtypes.push_back(origt);
    }
  }

  if (receptor_coords.size() == 0) { // Do once at setup

    // Reserve memory, but size() == 0
    receptor_coords.reserve(n);
    receptor_map.reserve(n);

    // Append flex
    const atomv &atoms = m.get_movable_atoms();
    const vecv &coords = m.coordinates();
    for (sz i = 0; i < num_flex_atoms; i++) {
      if (!atoms[i].iscov) {
        const vec &c = coords[i];
        receptor_coords.push_back(float3({c[0], c[1], c[2]}));
        receptor_map.push_back(i);
      }
    }

    // Append inflex
    for (sz i = m.m_num_movable_atoms, n = coords.size(); i < n; i++) {
      if (!atoms[i].iscov) {
        const vec &c = coords[i];
        receptor_coords.push_back(float3({c[0], c[1], c[2]}));
        receptor_map.push_back(i);
      }
    }

    // Append rigid receptor
    auto cbegin_rigid = m.get_fixed_atoms().cbegin();
    auto cend_rigid = m.get_fixed_atoms().cend();
    std::transform(cbegin_rigid, cend_rigid, std::back_inserter(receptor_coords), [](const atom &a) -> float3 {
      const vec &coord = a.coords;
      return float3({coord[0], coord[1], coord[2]});
    });

  } else if (receptor_coords.size() == n) { // Update flex coordinates at every call
    const atomv &atoms = m.get_movable_atoms();
    const vecv &coords = m.coordinates();
    for (sz i = 0; i < num_flex_atoms; i++) {
      if (!atoms[i].iscov) {
        const vec &c = coords[i];
        receptor_coords[i] = float3({c[0], c[1], c[2]});
      }
    }
  }
}

// reset center to be around ligand; reset receptor transformation
// call this before minimizing a ligand
void DLScorer::set_center_from_model(model &m) {

  // reset protein - not actually needed any more
  m.rec_conf.position = vec(0, 0, 0);
  m.rec_conf.orientation = qt(1, 0, 0, 0);

  // calc center for ligand
  current_center = vec(0, 0, 0);
  for (auto coord : m.coordinates()) {
    current_center += coord;
  }
  current_center /= (float)m.coordinates().size();

  if (cnnopts.verbose) {
    std::cout << "new center: ";
    current_center.print(std::cout);
    std::cout << "\n";
  }
}
