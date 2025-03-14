[![codecov](https://codecov.io/gh/gnina/gnina/branch/master/graph/badge.svg?token=kOC9GUKEGa)](https://codecov.io/gh/gnina/gnina)
[![Github CI](https://github.com/gnina/gnina/actions/workflows/CI.yml/badge.svg)](https://github.com/gnina/gnina/actions)

gnina (pronounced NEE-na) is a molecular docking program with integrated support for scoring and optimizing ligands using convolutional neural networks. It is a fork of [smina](http://smina.sf.net/), which is a fork of [AutoDock Vina](http://vina.scripps.edu/).  

Help
====
 
Please [subscribe to our slack team](https://join.slack.com/t/gninacnn/shared_invite/enQtNTY3ODk2ODk5OTU5LTkzMjY1ZTE3YjJlZmIxOWI2OTU3Y2RlMTIyYmM2YmFmYTU1NTk5ZTBmMjUwMGRhYzk1ZjY5N2E4Y2I5YWU5YWI).
An example colab notebook showing how to use gnina is available <a href="https://colab.research.google.com/drive/1QYo5QLUE80N_G28PlpYs6OKGddhhd931?usp=sharing">here</a>. We also hosted a workshp on using gnina ([video](https://www.youtube.com/watch?v=MG3Srzi5kZ0), [slides](https://bits.csb.pitt.edu/rsc_workshop2021/docking_with_gnina.slides.html#/)).

Citation
========
If you find gnina useful, please cite our paper(s):  

**GNINA 1.0: Molecular docking with deep learning** (Primary application citation)  
A McNutt, P Francoeur, R Aggarwal, T Masuda, R Meli, M Ragoza, J Sunseri, DR Koes. *J. Cheminformatics*, 2021  
[link](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2) [PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8191141/) [ChemRxiv](https://chemrxiv.org/articles/preprint/GNINA_1_0_Molecular_Docking_with_Deep_Learning/13578140)

**Protein–Ligand Scoring with Convolutional Neural Networks**  (Primary methods citation)  
M Ragoza, J Hochuli, E Idrobo, J Sunseri, DR Koes. *J. Chem. Inf. Model*, 2017  
[link](http://pubs.acs.org/doi/full/10.1021/acs.jcim.6b00740) [PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5479431/) [arXiv](https://arxiv.org/abs/1612.02751)  

**Ligand pose optimization with atomic grid-based convolutional neural networks**  
M Ragoza, L Turner, DR Koes. *Machine Learning for Molecules and Materials NIPS 2017 Workshop*, 2017  
[arXiv](https://arxiv.org/abs/1710.07400)  

**Visualizing convolutional neural network protein-ligand scoring**  
J Hochuli, A Helbling, T Skaist, M Ragoza, DR Koes.  *Journal of Molecular Graphics and Modelling*, 2018  
[link](https://www.sciencedirect.com/science/article/pii/S1093326318301670) [PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6343664/) [arXiv](https://arxiv.org/abs/1803.02398)

**Convolutional neural network scoring and minimization in the D3R 2017 community challenge**  
J Sunseri, JE King, PG Francoeur, DR Koes.  *Journal of computer-aided molecular design*, 2018  
[link](https://link.springer.com/article/10.1007/s10822-018-0133-y) [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/29992528)

**Three-Dimensional Convolutional Neural Networks and a Cross-Docked Data Set for Structure-Based Drug Design**  
PG Francoeur, T Masuda, J Sunseri, A Jia, RB Iovanisci, I Snyder, DR Koes. *J. Chem. Inf. Model*, 2020  
[link](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00411) [PubMed](https://pubmed.ncbi.nlm.nih.gov/32865404/) [Chemrxiv](https://chemrxiv.org/articles/preprint/3D_Convolutional_Neural_Networks_and_a_CrossDocked_Dataset_for_Structure-Based_Drug_Design/11833323/1)

**Virtual Screening with Gnina 1.0**
J Sunseri, DR Koes D. *Molecules*, 2021
[link](https://www.mdpi.com/1420-3049/26/23/7369) [Preprints](https://www.preprints.org/manuscript/202111.0329/v1)

Docker
======

A pre-built docker image is available [here](https://hub.docker.com/u/gnina) and Dockerfiles are [here](https://github.com/gnina/gnina/tree/master/docker).

Installation
============

We recommend that you use the [pre-built binary](https://github.com/gnina/gnina/releases/latest) unless you have significant experience building software on Linux, in which case building from source might result in an executable more optimized for your system.

### Ubuntu 22.04
```
apt-get  install build-essential git cmake wget libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest libjsoncpp-dev

```

[Follow NVIDIA's instructions](http://docs.nvidia.com/cuda/cuda-installation-guide-linux/#axzz4TWipdwX1) to install the latest version of CUDA (>= 12.0 is required). **Make sure `nvcc` is in your PATH.**

*Optionally* install [cuDNN](https://developer.nvidia.com/rdp/cudnn-archive).

#
Install OpenBabel3.  Note there are errors in bond order determination in version 3.1.1 and older.
```
git clone https://github.com/openbabel/openbabel.git
cd openbabel
mkdir build
cd build
cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON ..
make
make install
```

#
Install gnina
```
git clone https://github.com/gnina/gnina.git
cd gnina
mkdir build
cd build
cmake ..
make
make install
```

### [WSL2 Ubuntu 22.04](https://github.com/gnina/gnina/issues/247)
```bash
sudo apt-get remove nvidia-cuda-toolkit
wget https://developer.download.nvidia.com/compute/cuda/12.4.0/local_installers/cuda_12.4.0_550.54.14_linux.run
chmod 700 cuda_12.4.0_550.54.14_linux.run
sudo sh cuda_12.4.0_550.54.14_linux.run
wget https://developer.download.nvidia.com/compute/cudnn/9.0.0/local_installers/cudnn-local-repo-ubuntu2204-9.0.0_1.0-1_amd64.deb
sudo dpkg -i cudnn-local-repo-ubuntu2204-9.0.0_1.0-1_amd64.deb
sudo cp /var/cudnn-local-repo-ubuntu2204-9.0.0/cudnn-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cudnn-cuda-12
apt-get install build-essential git cmake wget libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest libjsoncpp-dev

git clone https://github.com/openbabel/openbabel.git
cd openbabel
mkdir build
cd build
cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON ..
make -j8
sudo make install

git clone https://github.com/gnina/gnina.git
cd gnina
mkdir build
cd build
cmake ..
make -j8
sudo make install
```

#
If you are building for systems with different GPUs (e.g. in a cluster environment), configure with `-DCMAKE_CUDA_ARCHITECTURES=all`.   
Note that the cmake build will automatically fetch and install [libmolgrid](https://github.com/gnina/libmolgrid) and [torch](https://github.com/pytorch/pytorch) if they are not already installed.


The scripts provided in `gnina/scripts` have additional python dependencies that must be installed. 

Usage
=====

To dock ligand `lig.sdf` to a binding site on `rec.pdb` defined by another ligand `orig.sdf`:
```
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf -o docked.sdf.gz
```

To perform docking with flexible sidechain residues within 3.5 Angstroms of `orig.sdf` (generally not recommend unless prior knowledge indicates pocket is highly flexible):
```
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf --flexdist_ligand orig.sdf --flexdist 3.5 -o flex_docked.sdf.gz
```

To perform whole protein docking:
```
gnina -r rec.pdb -l lig.sdf --autobox_ligand rec.pdb -o whole_docked.sdf.gz --exhaustiveness 64
```

To utilize the default ensemble CNN in the energy minimization during the refinement step of docking (10 times slower than the default rescore option):
```
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf --cnn_scoring refinement -o cnn_refined.sdf.gz
```

To utilize the default ensemble CNN for every step of docking (1000 times slower than the default rescore option):
```
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf --cnn_scoring all -o cnn_all.sdf.gz
```

To utilize all empirical scoring using the [Vinardo](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155183) scoring function:
```
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf --scoring vinardo --cnn_scoring none -o vinardo_docked.sdf.gz
```

To utilize a different CNN during docking (see help for possible options):
```

gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf --cnn dense -o dense_docked.sdf.gz
```

To minimize and score ligands `ligs.sdf` already positioned in a binding site:
```
gnina -r rec.pdb -l ligs.sdf --minimize -o minimized.sdf.gz
```

To covalently dock a pyrazole to a specific iron atom on the receptor with the bond formed between a nitrogen of the pyrazole and the iron.
```
gnina  -r rec.pdb.gz -l conformer.sdf.gz --autobox_ligand bindingsite.sdf.gz --covalent_rec_atom A:601:FE --covalent_lig_atom_pattern '[$(n1nccc1)]' -o output.sdf.gz 
```

The same as above, but with the covalently bonding ligand atom manually positioned (instead of using OpenBabel binding heuristics) and the ligand/residue complex UFF optimized.
```
gnina  -r rec.pdb.gz -l conformer.sdf.gz --autobox_ligand bindingsite.sdf.gz --covalent_lig_atom_position -11.796,31.887,72.682  --covalent_optimize_lig  --covalent_rec_atom A:601:FE --covalent_lig_atom_pattern '[$(n1nccc1)]' -o output.sdf.gz 
```

All options:
```
Input:
  -r [ --receptor ] arg              rigid part of the receptor
  --flex arg                         flexible side chains, if any (PDBQT)
  -l [ --ligand ] arg                ligand(s)
  --flexres arg                      flexible side chains specified by comma 
                                     separated list of chain:resid
  --flexdist_ligand arg              Ligand to use for flexdist
  --flexdist arg                     set all side chains within specified 
                                     distance to flexdist_ligand to flexible
  --flex_limit arg                   Hard limit for the number of flexible 
                                     residues
  --flex_max arg                     Retain at at most the closest flex_max 
                                     flexible residues

Search space (required):
  --center_x arg                     X coordinate of the center
  --center_y arg                     Y coordinate of the center
  --center_z arg                     Z coordinate of the center
  --size_x arg                       size in the X dimension (Angstroms)
  --size_y arg                       size in the Y dimension (Angstroms)
  --size_z arg                       size in the Z dimension (Angstroms)
  --autobox_ligand arg               Ligand to use for autobox
  --autobox_add arg                  Amount of buffer space to add to 
                                     auto-generated box (default +4 on all six 
                                     sides)
  --autobox_extend arg (=1)          Expand the autobox if needed to ensure the
                                     input conformation of the ligand being 
                                     docked can freely rotate within the box.
  --no_lig                           no ligand; for sampling/minimizing 
                                     flexible residues

Covalent docking:
  --covalent_rec_atom arg            Receptor atom ligand is covalently bound 
                                     to.  Can be specified as 
                                     chain:resnum:atom_name or as x,y,z 
                                     Cartesian coordinates.
  --covalent_lig_atom_pattern arg    SMARTS expression for ligand atom that 
                                     will covalently bind protein.
  --covalent_lig_atom_position arg   Optional.  Initial placement of covalently
                                     bonding ligand atom in x,y,z Cartesian 
                                     coordinates.  If not specified, 
                                     OpenBabel's GetNewBondVector function will
                                     be used to position ligand.
  --covalent_fix_lig_atom_position   If covalent_lig_atom_position is 
                                     specified, fix the ligand atom to this 
                                     position as opposed to using this position
                                     to define the initial structure.
  --covalent_bond_order arg (=1)     Bond order of covalent bond. Default 1.
  --covalent_optimize_lig            Optimize the covalent complex of ligand 
                                     and residue using UFF. This will change 
                                     bond angles and lengths of the ligand.

Scoring and minimization options:
  --scoring arg                      specify alternative built-in scoring 
                                     function: ad4_scoring default dkoes_fast 
                                     dkoes_scoring dkoes_scoring_old vina 
                                     vinardo
  --custom_scoring arg               custom scoring function file
  --custom_atoms arg                 custom atom type parameters file
  --score_only                       score provided ligand pose
  --local_only                       local search only using autobox (you 
                                     probably want to use --minimize)
  --minimize                         energy minimization
  --randomize_only                   generate random poses, attempting to avoid
                                     clashes
  --num_mc_steps arg                 fixed number of monte carlo steps to take 
                                     in each chain
  --max_mc_steps arg                 cap on number of monte carlo steps to take
                                     in each chain
  --num_mc_saved arg                 number of top poses saved in each monte 
                                     carlo chain
  --temperature arg                  temperature for metropolis accept 
                                     criterion
  --minimize_iters arg (=0)          number iterations of steepest descent; 
                                     default scales with rotors and usually 
                                     isn't sufficient for convergence
  --accurate_line                    use accurate line search
  --simple_ascent                    use simple gradient ascent
  --minimize_early_term              Stop minimization before convergence 
                                     conditions are fully met.
  --minimize_single_full             During docking perform a single full 
                                     minimization instead of a truncated 
                                     pre-evaluate followed by a full.
  --approximation arg                approximation (linear, spline, or exact) 
                                     to use
  --factor arg                       approximation factor: higher results in a 
                                     finer-grained approximation
  --force_cap arg                    max allowed force; lower values more 
                                     gently minimize clashing structures
  --user_grid arg                    Autodock map file for user grid data based
                                     calculations
  --user_grid_lambda arg (=-1)       Scales user_grid and functional scoring
  --print_terms                      Print all available terms with default 
                                     parameterizations
  --print_atom_types                 Print all available atom types

Convolutional neural net (CNN) scoring:
  --cnn_scoring arg (=1)             Amount of CNN scoring: none, rescore 
                                     (default), refinement, metrorescore 
                                     (metropolis+rescore), metrorefine 
                                     (metropolis+refine), all
  --cnn arg                          built-in model to use, specify 
                                     PREFIX_ensemble to evaluate an ensemble of
                                     models starting with PREFIX: 
                                     all_default_to_default_1_3_1 
                                     all_default_to_default_1_3_2 
                                     all_default_to_default_1_3_3 
                                     crossdock_default2018 
                                     crossdock_default2018_1 
                                     crossdock_default2018_1_3 
                                     crossdock_default2018_1_3_1 
                                     crossdock_default2018_1_3_2 
                                     crossdock_default2018_1_3_3 
                                     crossdock_default2018_1_3_4 
                                     crossdock_default2018_2 
                                     crossdock_default2018_3 
                                     crossdock_default2018_4 
                                     crossdock_default2018_KD_1 
                                     crossdock_default2018_KD_2 
                                     crossdock_default2018_KD_3 
                                     crossdock_default2018_KD_4 
                                     crossdock_default2018_KD_5 default1.0 
                                     default2017 dense dense_1 dense_1_3 
                                     dense_1_3_1 dense_1_3_2 dense_1_3_3 
                                     dense_1_3_4 dense_1_3_PT_KD 
                                     dense_1_3_PT_KD_1 dense_1_3_PT_KD_2 
                                     dense_1_3_PT_KD_3 dense_1_3_PT_KD_4 
                                     dense_1_3_PT_KD_def2018 
                                     dense_1_3_PT_KD_def2018_1 
                                     dense_1_3_PT_KD_def2018_2 
                                     dense_1_3_PT_KD_def2018_3 
                                     dense_1_3_PT_KD_def2018_4 dense_2 dense_3 
                                     dense_4 fast general_default2018 
                                     general_default2018_1 
                                     general_default2018_2 
                                     general_default2018_3 
                                     general_default2018_4 
                                     general_default2018_KD_1 
                                     general_default2018_KD_2 
                                     general_default2018_KD_3 
                                     general_default2018_KD_4 
                                     general_default2018_KD_5 
                                     redock_default2018 redock_default2018_1 
                                     redock_default2018_1_3 
                                     redock_default2018_1_3_1 
                                     redock_default2018_1_3_2 
                                     redock_default2018_1_3_3 
                                     redock_default2018_1_3_4 
                                     redock_default2018_2 redock_default2018_3 
                                     redock_default2018_4 redock_default2018_KD
                                     _1 redock_default2018_KD_2 
                                     redock_default2018_KD_3 
                                     redock_default2018_KD_4 
                                     redock_default2018_KD_5
  --cnn_model arg                    torch cnn model file; if not specified a 
                                     default model ensemble will be used
  --cnn_rotation arg (=0)            evaluate multiple rotations of pose (max 
                                     24)
  --cnn_mix_emp_force                Merge CNN and empirical minus forces
  --cnn_mix_emp_energy               Merge CNN and empirical energy
  --cnn_empirical_weight arg (=1)    Weight for scaling and merging empirical 
                                     force and energy 
  --cnn_outputdx                     Dump .dx files of atom grid gradient.
  --cnn_outputxyz                    Dump .xyz files of atom gradient.
  --cnn_xyzprefix arg (=gradient)    Prefix for atom gradient .xyz files
  --cnn_center_x arg                 X coordinate of the CNN center
  --cnn_center_y arg                 Y coordinate of the CNN center
  --cnn_center_z arg                 Z coordinate of the CNN center
  --cnn_verbose                      Enable verbose output for CNN debugging

Output:
  -o [ --out ] arg                   output file name, format taken from file 
                                     extension
  --out_flex arg                     output file for flexible receptor residues
  --log arg                          optionally, write log file
  --atom_terms arg                   optionally write per-atom interaction term
                                     values
  --atom_term_data                   embedded per-atom interaction terms in 
                                     output sd data
  --pose_sort_order arg (=0)         How to sort docking results: CNNscore 
                                     (default), CNNaffinity, Energy
  --full_flex_output                 Output entire structure for out_flex, not 
                                     just flexible residues.

Misc (optional):
  --cpu arg                          the number of CPUs to use (the default is 
                                     to try to detect the number of CPUs or, 
                                     failing that, use 1)
  --seed arg                         explicit random seed
  --exhaustiveness arg (=8)          exhaustiveness of the global search 
                                     (roughly proportional to time)
  --num_modes arg (=9)               maximum number of binding modes to 
                                     generate
  --min_rmsd_filter arg (=1)         rmsd value used to filter final poses to 
                                     remove redundancy
  -q [ --quiet ]                     Suppress output messages
  --addH arg                         automatically add hydrogens in ligands (on
                                     by default)
  --stripH arg                       remove hydrogens from molecule _after_ 
                                     performing atom typing for efficiency (off
                                     by default)
  --device arg (=0)                  GPU device to use
  --no_gpu                           Disable GPU acceleration, even if 
                                     available.

Configuration file (optional):
  --config arg                       the above options can be put here

Information (optional):
  --help                             display usage summary
  --help_hidden                      display usage summary with hidden options
  --version                          display program version
```

CNN Scoring
===========

`--cnn_scoring` determines at what points of the docking procedure that the CNN scoring function is used.
 * `none` - No CNNs used for docking. Uses the specified empirical scoring function throughout.
 * `rescore` (default) - CNN used for reranking of final poses. Least computationally expensive CNN option.
 * `refinement` - CNN used to refine poses after Monte Carlo chains and for final ranking of output poses. 10x slower than `rescore` when using a GPU.
 * `all` - CNN used as the scoring function throughout the whole procedure. Extremely computationally intensive and not recommended.

The default CNN scoring function is an ensemble of 3 models selected to balance pose prediction performance and runtime: dense_1_3, dense_1_3_PT_KD_3, crossdock_default2018_KD_4.  The GNINA 1.0 default ensemble is still available as `default1.0`.   More information on these various models can be found in the papers listed above.

Training
========

Scripts for training pytorch GNINA models and pretrained models can found at [https://github.com/RMeli/gnina-torch](https://github.com/RMeli/gnina-torch).
Example code for converting a pytorch model into a gnina usable model file is shown below.  The metadata should provide information about the input grid resolution, dimension and atom typing.  If not provided, defaults will be used.

```python
    d = {
        'resolution': 0.5,
        'dimension' : 23.5,
        'recmap' : '''AliphaticCarbonXSHydrophobe 
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
    GenericMetal Boron Manganese Magnesium Iron''',
        
    'ligmap': '''AliphaticCarbonXSHydrophobe 
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
    GenericMetal Boron Manganese Magnesium Zinc Calcium Iron'''
    }
    
    extra = {'metadata':json.dumps(d)}
    z = torch.zeros((1,28,48,48,48))

    script = torch.jit.trace(model, z)
    script.save('gnina_model.pt',_extra_files=extra)    
```

Legacy scripts for training Caffe models can be found at [https://github.com/gnina/scripts](https://github.com/gnina/scripts)
and sample models at [https://github.com/gnina/models](https://github.com/gnina/models).


The DUD-E docked poses used in the original paper can be found [here](http://bits.csb.pitt.edu/files/docked_dude.tar), but we [do not recommend](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0220113) training virtual screening models on DUD-E. The CrossDocked2020 set is [here](https://github.com/gnina/models/tree/master/data/CrossDocked2020).

License
=======

gnina is dual licensed under GPL and Apache.  The GPL license is necessitated by the use of OpenBabel (which is GPL licensed).  In order to use gnina under the Apache license only, all references to OpenBabel must be removed from the source code.
