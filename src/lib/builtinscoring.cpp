#include "builtinscoring.h"

builtin_scoring::builtin_scoring()
{
	//set all builtin functions

	add("vina","gauss(o=0,_w=0.5,_c=8)", -0.035579);
	add("vina","gauss(o=3,_w=2,_c=8)", -0.005156);
	add("vina","repulsion(o=0,_c=8)", 0.840245);
	add("vina","hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069);
	add("vina","non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439);
	add("vina","num_tors_div", 5 * 0.05846 / 0.1 - 1);

	functions["default"] = functions["vina"];

	add("dkoes_scoring","vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.009900);
	add("dkoes_scoring","non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.153055);
	add("dkoes_scoring","ad4_solvation(d-sigma=3.6,_s/q=0.01097,_c=8)", 0.048934);
	add("dkoes_scoring","num_tors_sqr", 0.317267);
	add("dkoes_scoring","constant_term", -2.469020);
	/* trained with openbabel partial charges
	 weights.push_back(0.010764); //vdw
	 weights.push_back(-0.156861); //hbond
	 weights.push_back(0.062407); //desolvation
	 weights.push_back(0.337036); //tors sqr
	 weights.push_back(-2.231827); //constant
	 */

	add("dkoes_scoring_old","vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.010607);
	add("dkoes_scoring_old","non_dir_h_bond(g=-0.7,_b=0,_c=8)", 0.197201);
	add("dkoes_scoring_old","num_tors_sqr", .285035);
	add("dkoes_scoring_old","constant_term", -2.585651);

	add("dkoes_fast","vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.008962);
	add("dkoes_fast","non_dir_h_bond(g=-0.7,_b=0,_c=8)", 0.387739);
	add("dkoes_fast","num_tors_sqr", .285035);
	add("dkoes_fast","constant_term", -2.467357);

	add("ad4_scoring","vdw(i=6,_j=12,_s=0,_^=100,_c=8)", 0.1560);
	add("ad4_scoring","non_dir_h_bond_lj(o=-0.7,_^=100,_c=8)", -0.0974);
	add("ad4_scoring","ad4_solvation(d-sigma=3.5,_s/q=0.01097,_c=8)", 0.1159);
	add("ad4_scoring","electrostatic(i=1,_^=100,_c=8)", 0.1465);
	add("ad4_scoring","num_tors_add", 0.2744);
}

void builtin_scoring::print_functions(std::ostream& out)
{
	std::vector<std::string> names;
	for(funcmap::iterator itr = functions.begin(), en = functions.end(); itr != en; ++itr)
	{
		names.push_back(itr->first);
	}
	std::sort(names.begin(),names.end());

	for(unsigned i = 0, n = names.size(); i < n; i++)
	{
		out << names[i] << "\n";
	}
}

//global object for getting scoring functions
builtin_scoring builtin_scoring_functions;
