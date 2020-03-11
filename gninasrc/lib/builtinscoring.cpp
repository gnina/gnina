#include "builtinscoring.h"


namespace smina_atom_type
{
//vinardo scoring has its own atom parameters
const info vinardo_data[NumTypes] = { //el, ad, xs
    {Hydrogen, "Hydrogen", "H", 1,  1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.000000, false,  false,  false,  false},
    {PolarHydrogen, "PolarHydrogen", "HD", 1, 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.000000, false,  false,  false,  false},
    {AliphaticCarbonXSHydrophobe, "AliphaticCarbonXSHydrophobe", "C", 6,  2.000000, 0.150000, -0.001430,  33.510300,  0.770000, 2.000000, true, false,  false,  false},
    {AliphaticCarbonXSNonHydrophobe, "AliphaticCarbonXSNonHydrophobe", "C", 6, 2.000000, 0.150000, -0.001430,  33.510300,  0.770000, 2.000000, false,  false,  false,  false},
    {AromaticCarbonXSHydrophobe, "AromaticCarbonXSHydrophobe", "A", 6, 2.000000, 0.150000, -0.000520,  33.510300,  0.770000, 1.900000, true, false,  false,  false},
    {AromaticCarbonXSNonHydrophobe, "AromaticCarbonXSNonHydrophobe", "A", 6, 2.000000, 0.150000, -0.000520,  33.510300,  0.770000, 1.900000, true,  false,  false,  false},
    {Nitrogen, "Nitrogen", "N", 7,  1.750000, 0.160000, -0.001620,  22.449300,  0.750000, 1.700000, false,  false,  false,  true},
    {NitrogenXSDonor, "NitrogenXSDonor", "N", 7, 1.750000, 0.160000, -0.001620,  22.449300,  0.750000, 1.700000, false,  true, false,  true},
    {NitrogenXSDonorAcceptor,"NitrogenXSDonorAcceptor", "NA", 7, 1.750000, 0.160000, -0.001620,  22.449300,  0.750000, 1.700000, false,  true, true, true},
    {NitrogenXSAcceptor, "NitrogenXSAcceptor", "NA", 7, 1.750000, 0.160000, -0.001620,  22.449300,  0.750000, 1.700000, false,  false,  true, true},
    {Oxygen, "Oxygen", "O", 8, 1.600000, 0.200000, -0.002510,  17.157300,  0.730000, 1.600000, false,  false,  false,  true},
    {OxygenXSDonor,"OxygenXSDonor", "O", 8, 1.600000, 0.200000, -0.002510,  17.157300,  0.730000, 1.600000, false,  true, false,  true},
    {OxygenXSDonorAcceptor,"OxygenXSDonorAcceptor", "OA", 8, 1.600000, 0.200000, -0.002510,  17.157300,  0.730000, 1.600000, false,  true, true, true},
    {OxygenXSAcceptor, "OxygenXSAcceptor", "OA", 8, 1.600000, 0.200000, -0.002510,  17.157300,  0.730000, 1.600000, false,  false,  true, true},
    {Sulfur,"Sulfur", "S", 16, 2.000000, 0.200000, -0.002140,  33.510300,  1.020000, 2.000000, false,  false,  false,  true},
    {SulfurAcceptor,"SulfurAcceptor", "SA", 16, 2.000000, 0.200000, -0.002140,  33.510300,  1.020000, 2.000000, true,  false,  false,  true},
    {Phosphorus, "Phosphorus", "P", 15, 2.100000, 0.200000, -0.001100,  38.792400,  1.060000, 2.100000, false,  false,  false,  true},
    {Fluorine, "Fluorine", "F", 9, 1.545000, 0.080000, -0.001100,  15.448000,  0.710000, 1.500000, true, false,  false,  true},
    {Chlorine,"Chlorine", "Cl", 17, 2.045000, 0.276000, -0.001100,  35.823500,  0.990000, 1.800000, true, false,  false,  true},
    {Bromine,"Bromine", "Br", 35, 2.165000, 0.389000, -0.001100,  42.566100,  1.140000, 2.000000, true, false,  false,  true},
    {Iodine,"Iodine", "I", 53, 2.360000, 0.550000, -0.001100,  55.058500,  1.330000, 2.200000, true, false,  false,  true},
    {Magnesium,"Magnesium", "Mg", 12, 0.650000, 0.875000, -0.001100,  1.560000, 1.300000, 1.200000, false,  true, false,  true},
    {Manganese,"Manganese", "Mn", 25, 0.650000, 0.875000, -0.001100,  2.140000, 1.390000, 1.200000, false,  true, false,  true},
    {Zinc, "Zinc","Zn", 30, 0.740000, 0.550000, -0.001100,  1.700000, 1.310000, 1.200000, false,  true, false,  true},
    {Calcium,"Calcium", "Ca", 20, 0.990000, 0.550000, -0.001100,  2.770000, 1.740000, 1.200000, false,  true, false,  true},
    {Iron,"Iron", "Fe", 26, 0.650000, 0.010000, -0.001100,  1.840000, 1.250000, 1.200000, false,  true, false,  true},
    {GenericMetal,"GenericMetal", "M", 0, 1.200000, 0.000000, -0.001100,  22.449300,  1.750000, 1.200000, false,  true, false,  true}
};
} //smina_atom_type namespace

builtin_scoring::builtin_scoring() {
  //set all builtin functions

  add("vina", "gauss(o=0,_w=0.5,_c=8)", -0.035579);
  add("vina", "gauss(o=3,_w=2,_c=8)", -0.005156);
  add("vina", "repulsion(o=0,_c=8)", 0.840245);
  add("vina", "hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069);
  add("vina", "non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439);
  add("vina", "num_tors_div", 5 * 0.05846 / 0.1 - 1);

  functions["default"] = functions["vina"];

  add("vinardo","gauss(o=0,_w=0.8,_c=8)", -0.045);
  add("vinardo","repulsion(o=0,_c=8)", 0.80);
  add("vinardo","hydrophobic(g=0.0,_b=2.5,_c=8)", -0.035);
  add("vinardo","non_dir_h_bond(g=-0.6,_b=0,_c=8)", -0.60);
  add("vinardo","num_tors_div", 5 * 0.02 / 0.1 - 1);
  addparams("vinardo", smina_atom_type::vinardo_data);

  add("dkoes_scoring", "vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.009900);
  add("dkoes_scoring", "non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.153055);
  add("dkoes_scoring", "ad4_solvation(d-sigma=3.6,_s/q=0.01097,_c=8)",
      0.048934);
  add("dkoes_scoring", "num_tors_sqr", 0.317267);
  add("dkoes_scoring", "constant_term", -2.469020);
  /* trained with openbabel partial charges
   weights.push_back(0.010764); //vdw
   weights.push_back(-0.156861); //hbond
   weights.push_back(0.062407); //desolvation
   weights.push_back(0.337036); //tors sqr
   weights.push_back(-2.231827); //constant
   */

  add("dkoes_scoring_old", "vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.010607);
  add("dkoes_scoring_old", "non_dir_h_bond(g=-0.7,_b=0,_c=8)", 0.197201);
  add("dkoes_scoring_old", "num_tors_sqr", .285035);
  add("dkoes_scoring_old", "constant_term", -2.585651);

  add("dkoes_fast", "vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.008962);
  add("dkoes_fast", "non_dir_h_bond(g=-0.7,_b=0,_c=8)", 0.387739);
  add("dkoes_fast", "num_tors_sqr", .285035);
  add("dkoes_fast", "constant_term", -2.467357);

  add("ad4_scoring", "vdw(i=6,_j=12,_s=0,_^=100,_c=8)", 0.1560);
  add("ad4_scoring", "non_dir_h_bond_lj(o=-0.7,_^=100,_c=8)", 0.0974);
  add("ad4_scoring", "ad4_solvation(d-sigma=3.5,_s/q=0.01097,_c=8)", 0.1159);
  add("ad4_scoring", "electrostatic(i=1,_^=100,_c=8)", 0.1465);
  add("ad4_scoring", "num_tors_add", 0.2744);
}

void builtin_scoring::print_functions(std::ostream& out) {
  std::vector < std::string > names;
  for (funcmap::iterator itr = functions.begin(), en = functions.end();
      itr != en; ++itr) {
    names.push_back(itr->first);
  }
  std::sort(names.begin(), names.end());

  for (unsigned i = 0, n = names.size(); i < n; i++) {
    out << names[i] << "\n";
  }
}

//global object for getting scoring functions
builtin_scoring builtin_scoring_functions;
