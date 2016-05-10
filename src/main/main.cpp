#include <boost/program_options.hpp>

#include <iostream>
#include <string>
#include <exception>
#include <vector> // ligand paths#include <cmath> // for ceila
#include <algorithm>
#include <iterator>#include <boost/filesystem/fstream.hpp>#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?#include <boost/lexical_cast.hpp>
#include <boost/assign.hpp>
#include "parse_pdbqt.h"#include "parallel_mc.h"
#include "file.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "non_cache_gpu.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "quasi_newton.h"
#include "tee.h"
#include "custom_terms.h"
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
#include <boost/shared_ptr.hpp>
#include "coords.h"
#include "obmolopener.h"
#include "gpucode.h"
#include "precalculate_gpu.h"
#include <boost/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>
#include "array3d.h"
#include "grid.h"
#include "molgetter.h"
#include "result_info.h"
#include "box.h"
#include "flexinfo.h"
#include "builtinscoring.h"

#include <cuda_profiler_api.h>

using namespace boost::iostreams;
using boost::filesystem::path;


//just a collection of user-specified configurations
struct user_settings
{
	fl energy_range;
	sz num_modes;
	fl out_min_rmsd;
	fl forcecap;
	int seed;
	int verbosity;
	int cpu;
	int exhaustiveness;
	bool score_only;
	bool randomize_only;
	bool local_only;
	bool dominimize;
	bool include_atom_info;

	//reasonable defaults
	user_settings(): energy_range(2.0), num_modes(9), out_min_rmsd(1),
			forcecap(1000),seed(auto_seed()),verbosity(1), cpu(1), exhaustiveness(10),
			score_only(false), randomize_only(false), local_only(false),
			dominimize(false), include_atom_info(false)
	{

	}
};

void doing(int verbosity, const std::string& str, tee& log)
{
	if (verbosity > 1)
	{
		log << str << std::string(" ... ");
		log.flush();
	}
}

void done(int verbosity, tee& log)
{
	if (verbosity > 1)
	{
		log << "done.";
		log.endl();
	}
}
std::string default_output(const std::string& input_name)
{
	std::string tmp = input_name;
	if (tmp.size() >= 6 && tmp.substr(tmp.size() - 6, 6) == ".pdbqt")
		tmp.resize(tmp.size() - 6); // FIXME?
	return tmp + "_out.pdbqt";
}

void write_all_output(model& m, const output_container& out, sz how_many,
		std::ostream& outstream)
{
	if (out.size() < how_many)
		how_many = out.size();
	VINA_FOR(i, how_many)
	{
		m.set(out[i].c);
		m.write_model(outstream, i + 1); // so that model numbers start with 1
	}
}

//set m to a random conformer
fl do_randomization(model& m, const vec& corner1,
		const vec& corner2, int seed, int verbosity, tee& log)
{
	conf init_conf = m.get_initial_conf();
	rng generator(static_cast<rng::result_type>(seed));
	if (verbosity > 1)
	{
		log << "Using random seed: " << seed;
		log.endl();
	}
	const sz attempts = 10000;
	conf best_conf = init_conf;
	fl best_clash_penalty = 0;
	VINA_FOR(i, attempts)
	{
		conf c = init_conf;
		c.randomize(corner1, corner2, generator);
		m.set(c);
		fl penalty = m.clash_penalty();
		if (i == 0 || penalty < best_clash_penalty)
		{
			best_conf = c;
			best_clash_penalty = penalty;
		}
	}
	m.set(best_conf);
	if (verbosity > 1)
	{
		log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
		log.endl();
	}
	return best_clash_penalty;
}

void refine_structure(model& m, const precalculate& prec, non_cache& nc,
		output_type& out, const vec& cap, const minimization_params& minparm,
		grid& user_grid)
{
	change g(m.get_size());

	quasi_newton quasi_newton_par(minparm);
	const fl slope_orig = nc.getSlope();
	//try 5 times to get ligand into box
	//dkoes - you don't need a very strong constraint to keep ligands in the box,
	//but this factor can really bias the energy landscape
	fl slope = 10;
	VINA_FOR(p, 5)
	{
		nc.setSlope(slope);
		quasi_newton_par(m, prec, nc, out, g, cap, user_grid); //quasi_newton operator
		m.set(out.c); // just to be sure
		if (nc.within(m)) {
			break;
		}
		slope *= 10;
	}
	out.coords = m.get_heavy_atom_movable_coords();
	if (!nc.within(m))
		out.e = max_fl;
	nc.setSlope(slope_orig);
}

std::string vina_remark(fl e, fl lb, fl ub)
{
	std::ostringstream remark;
	remark.setf(std::ios::fixed, std::ios::floatfield);
	remark.setf(std::ios::showpoint);
	remark << "REMARK VINA RESULT: " << std::setw(9) << std::setprecision(1)
			<< e << "  " << std::setw(9) << std::setprecision(3) << lb << "  "
			<< std::setw(9) << std::setprecision(3) << ub << '\n';
	return remark.str();
}

output_container remove_redundant(const output_container& in, fl min_rmsd)
{
	output_container tmp;
	VINA_FOR_IN(i, in)
		add_to_output_container(tmp, in[i], min_rmsd, in.size());
	return tmp;
}

//dkoes - return all energies and rmsds to original conf with result
void do_search(model& m, const boost::optional<model>& ref,
		const weighted_terms& sf, const precalculate& prec, const igrid& ig,
		non_cache& nc, // nc.slope is changed
		const vec& corner1, const vec& corner2,
		const parallel_mc& par, const user_settings& settings,
		bool compute_atominfo, tee& log,
		const terms *t, grid& user_grid, std::vector<result_info>& results)
{
	boost::timer time;

	precalculate_exact exact_prec(sf); //use exact computations for final score
	conf_size s = m.get_size();
	conf c = m.get_initial_conf();
	fl e = max_fl;
	fl rmsd = 0;
	const vec authentic_v(settings.forcecap, settings.forcecap, settings.forcecap); //small cap restricts initial movement from clash
	if (settings.score_only)
	{
		fl intramolecular_energy = m.eval_intramolecular(exact_prec,
				authentic_v, c);
		naive_non_cache nnc(&exact_prec); // for out of grid issues
		e = m.eval_adjusted(sf, exact_prec, nnc, authentic_v, c,
				intramolecular_energy, user_grid);

		log << "Affinity: " << std::fixed << std::setprecision(5) << e
				<< " (kcal/mol)";
		log.endl();
		std::vector<flv> atominfo;
		flv term_values = t->evale_robust(m);
		log << "Intramolecular energy: " << std::fixed << std::setprecision(5)
				<< intramolecular_energy << "\n";

		log
		<< "Term values, before weighting:\n";
		log << std::setprecision(5);
		log << "## " << boost::replace_all_copy(m.get_name()," ","_");

		VINA_FOR_IN(i, term_values)
		{
			log << ' ' << term_values[i];
		}

		conf_independent_inputs in(m);
		const flv nonweight(1, 1.0);
		for (unsigned i = 0, n = t->conf_independent_terms.size(); i < n; i++)
		{
			flv::const_iterator pos = nonweight.begin();
			log << " " <<
					t->conf_independent_terms[i].eval(in, (fl) 0.0, pos);
		}
		log << '\n';

		results.push_back(result_info(e, -1, m));
		if (compute_atominfo)
			results.back().setAtomValues(m, &sf);
	}
	else if (settings.local_only)
	{
		vecv origcoords = m.get_heavy_atom_movable_coords();
		output_type out(c, e);
		doing(settings.verbosity, "Performing local search", log);
		refine_structure(m, prec, nc, out, authentic_v, par.mc.ssd_par.minparm,
				user_grid);
		done(settings.verbosity, log);

		//be as exact as possible for final score
		naive_non_cache nnc(&exact_prec); // for out of grid issues

		fl intramolecular_energy = m.eval_intramolecular(exact_prec,
				authentic_v,
				out.c);
		e = m.eval_adjusted(sf, exact_prec, nnc, authentic_v, out.c,
				intramolecular_energy, user_grid);

		vecv newcoords = m.get_heavy_atom_movable_coords();
		assert(newcoords.size() == origcoords.size());
		for (unsigned i = 0, n = newcoords.size(); i < n; i++)
		{
			rmsd += (newcoords[i] - origcoords[i]).norm_sqr();
		}
		rmsd /= newcoords.size();
		rmsd = sqrt(rmsd);
		log << "Affinity: " << std::fixed << std::setprecision(5) << e << "  "
				<< intramolecular_energy
				<< " (kcal/mol)\nRMSD: " << rmsd;
		;
		log.endl();
		if (!nc.within(m))
			log
			<< "WARNING: not all movable atoms are within the search space\n";

		m.set(out.c);
		done(settings.verbosity, log);
		results.push_back(result_info(e, rmsd, m));
		if (compute_atominfo)
			results.back().setAtomValues(m, &sf);
	}
	else
	{
		rng generator(static_cast<rng::result_type>(settings.seed));
		log << "Using random seed: " << settings.seed;
		log.endl();
		output_container out_cont;
		doing(settings.verbosity, "Performing search", log);
		par(m, out_cont, prec, ig, corner1, corner2, generator, user_grid);
		done(settings.verbosity, log);
		doing(settings.verbosity, "Refining results", log);
		VINA_FOR_IN(i, out_cont)
			refine_structure(m, prec, nc, out_cont[i], authentic_v,
					par.mc.ssd_par.minparm, user_grid);

		if (!out_cont.empty())
		{
			out_cont.sort();
			const fl best_mode_intramolecular_energy = m.eval_intramolecular(
					prec, authentic_v, out_cont[0].c);

			VINA_FOR_IN(i, out_cont)
				if (not_max(out_cont[i].e))
					out_cont[i].e = m.eval_adjusted(sf, prec, nc, authentic_v,
							out_cont[i].c, best_mode_intramolecular_energy,
							user_grid);
			// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
			out_cont.sort();
		}
		out_cont = remove_redundant(out_cont, settings.out_min_rmsd);

		done(settings.verbosity, log);

		log.setf(std::ios::fixed, std::ios::floatfield);
		log.setf(std::ios::showpoint);
		log << '\n';
		log << "mode |   affinity | dist from best mode\n";
		log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
		log << "-----+------------+----------+----------\n";

		model best_mode_model = m;
		if (!out_cont.empty())
			best_mode_model.set(out_cont.front().c);

		sz how_many = 0;
		VINA_FOR_IN(i, out_cont)
		{
			if (how_many >= settings.num_modes || !not_max(out_cont[i].e)
					|| out_cont[i].e > out_cont[0].e + settings.energy_range)
				break; // check energy_range sanity FIXME
			++how_many;
			log << std::setw(4) << i + 1 << "    " << std::setw(9)
					<< std::setprecision(1) << out_cont[i].e; // intermolecular_energies[i];
			m.set(out_cont[i].c);
			const model& r = ref ? ref.get() : best_mode_model;
			const fl lb = m.rmsd_lower_bound(r);
			const fl ub = m.rmsd_upper_bound(r);
			log << "  " << std::setw(9) << std::setprecision(3) << lb << "  "
					<< std::setw(9) << std::setprecision(3) << ub; // FIXME need user-readable error messages in case of failures

			log.endl();

			//dkoes - setup result_info
			results.push_back(result_info(out_cont[i].e, -1, m));
			if (compute_atominfo)
				results.back().setAtomValues(m, &sf);

		}
		done(settings.verbosity, log);

		if (how_many < 1)
		{
			log
					<< "WARNING: Could not find any conformations completely within the search space.\n"
					<< "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
			log.endl();
		}
	}
	//std::cout << "Refine time " << time.elapsed() << "\n";
}

void load_ent_values(const grid_dims& gd, std::istream& user_in,
		array3d<fl>& user_data)
{
	std::string line;
	user_data = array3d<fl>(gd[0].n + 1, gd[1].n + 1, gd[2].n + 1);

	for (sz z = 0; z < gd[2].n + 1; z++)
	{
		for (sz y = 0; y < gd[1].n + 1; y++)
		{
			for (sz x = 0; x < gd[0].n + 1; x++)
			{
				std::getline(user_in, line);
				user_data(x, y, z) = ::atof(line.c_str());
			}
		}
	}
	std::cout << user_data(gd[0].n - 3, gd[1].n, gd[2].n) << "\n";
}

void main_procedure(model& m, precalculate& prec,
		const boost::optional<model>& ref, // m is non-const (FIXME?)
		const user_settings& settings,
		bool no_cache, bool compute_atominfo, bool gpu_on,
		const grid_dims& gd, minimization_params minparm,
		const weighted_terms& wt, tee& log,
		std::vector<result_info>& results, grid& user_grid)
{
	doing(settings.verbosity, "Setting up the scoring function", log);

	done(settings.verbosity, log);

	vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
	vec corner2(gd[0].end, gd[1].end, gd[2].end);

	parallel_mc par;
	sz heuristic = m.num_movable_atoms()
			+ 10 * m.get_size().num_degrees_of_freedom();
	par.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
	par.mc.ssd_par.evals = unsigned((25 + m.num_movable_atoms()) / 3);
	if (minparm.maxiters == 0)
		minparm.maxiters = par.mc.ssd_par.evals;
	par.mc.ssd_par.minparm = minparm;
	par.mc.min_rmsd = 1.0;
	par.mc.num_saved_mins = settings.num_modes > 20 ? settings.num_modes : 20; //dkoes, support more than 20
	par.mc.hunt_cap = vec(10, 10, 10);
	par.num_tasks = settings.exhaustiveness;
	par.num_threads = settings.cpu;
	par.display_progress = true;

	szv_grid_cache gridcache(m, prec.cutoff_sqr());
	const fl slope = 1e6; // FIXME: too large? used to be 100
	if (settings.randomize_only)
	{
		fl e = do_randomization(m, corner1, corner2, settings.seed, settings.verbosity, log);
		results.push_back(result_info(e, -1, m));
		return;
	}
	else
	{

		non_cache *nc = NULL;
		if(gpu_on)
		{
			precalculate_gpu *gprec = dynamic_cast<precalculate_gpu*>(&prec);
			if(!gprec) abort();
			nc = new non_cache_gpu(gridcache, gd, gprec, slope);
		}
		else
		{
			nc = new non_cache(gridcache, gd, &prec, slope);
		}
    /* cudaProfilerStart(); */

		if (no_cache)
		{
			do_search(m, ref, wt, prec, *nc, *nc, corner1, corner2, par,
					settings, compute_atominfo, log,
					wt.unweighted_terms(), user_grid,
					results);
		}
		else
		{
			bool cache_needed = !(settings.score_only || settings.randomize_only || settings.local_only);
			if (cache_needed)
				doing(settings.verbosity, "Analyzing the binding site", log);
			cache c("scoring_function_version001", gd, slope);
			if (cache_needed)
			{
				std::vector<smt> atom_types_needed;
				m.get_movable_atom_types(atom_types_needed);
				c.populate(m, prec, atom_types_needed, user_grid);
			}
			if (cache_needed)
				done(settings.verbosity, log);
			do_search(m, ref, wt, prec, c, *nc, corner1, corner2, par,
					settings, compute_atominfo, log,
					wt.unweighted_terms(), user_grid, results);
		}
    /* cudaProfilerStop(); */

		delete nc;
	}
}


struct options_occurrence
{
	bool some;
	bool all;
	options_occurrence() :
			some(false), all(true)
	{
	} // convenience
	options_occurrence& operator+=(const options_occurrence& x)
	{
		some = some || x.some;
		all = all && x.all;
		return *this;
	}
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm,
		boost::program_options::options_description& d)
{
	options_occurrence tmp;
	VINA_FOR_IN(i, d.options())
	{
		const std::string& str = (*d.options()[i]).long_name();
		if ((str.substr(0, 4) == "size" || str.substr(0, 6) == "center"))
		{
			if (vm.count(str))
				tmp.some = true;
			else
				tmp.all = false;
		}
	}
	return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm,
		boost::program_options::options_description& d)
{
	VINA_FOR_IN(i, d.options())
	{
		const std::string& str = (*d.options()[i]).long_name();
		if ((str.substr(0, 4) == "size" || str.substr(0, 6) == "center")
				&& !vm.count(str))
			std::cerr << "Required parameter --" << str << " is missing!\n";
	}
}


template <class T>
inline void read_atomconstants_field(smina_atom_type::info& info, T (smina_atom_type::info::* field), unsigned int line, const std::string& field_name, std::istream& in)
{
	if (!(in >> (info.*field) ))
	{
		throw usage_error("Error at line " + boost::lexical_cast<std::string>(line) + " while reading field '" + field_name + "' from the atom constants file.");
	}
}

void setup_atomconstants_from_file(const std::string& atomconstants_file)
{
	std::ifstream file(atomconstants_file.c_str());
	if (file)
	{
		// create map from atom type names to indices
		boost::unordered_map<std::string, unsigned> atomindex;
		for(size_t i = 0u; i < smina_atom_type::NumTypes; ++i)
		{
			atomindex[smina_atom_type::default_data[i].smina_name] = i;
		}

		//parse each line of the file
		std::string line;
		unsigned lineno = 1;
		while(std::getline(file, line))
		{
			std::string name;
			std::stringstream ss(line);

			if(line.length() == 0 || line[0] == '#')
				continue;

			ss >> name;
			
			if(atomindex.count(name))
			{
				unsigned i = atomindex[name];
				smina_atom_type::info& info = smina_atom_type::data[i];

				//change this atom's parameters
#				define read_field(field) read_atomconstants_field(info, &smina_atom_type::info::field, lineno, #field, ss)
				read_field(ad_radius);
				read_field(ad_depth);
				read_field(ad_solvation);
				read_field(ad_volume);
				read_field(covalent_radius);
				read_field(xs_radius);
				read_field(xs_hydrophobe);
				read_field(xs_donor);
				read_field(xs_acceptor);
				read_field(ad_heteroatom);
#				undef read_field
			}
			else
			{
				std::cerr << "Line " << lineno << ": ommitting atom type name " << name << "\n";
			}
			lineno++;
		}
	}
	else
		throw usage_error("Error opening atom constants file:  " + atomconstants_file);
}

void print_atom_info(std::ostream& out)
{
	out << "#Name radius depth solvation volume covalent_radius xs_radius xs_hydrophobe xs_donor xs_acceptr ad_heteroatom\n";
	VINA_FOR(i, smina_atom_type::NumTypes) {
		smina_atom_type::info& info = smina_atom_type::data[i];
		out << info.smina_name;
		out << " " << info.ad_radius;
		out << " " << info.ad_depth;
		out << " " << info.ad_solvation;
		out << " " << info.ad_volume;
		out << " " << info.covalent_radius;
		out << " " << info.xs_radius;
		out << " " << info.xs_hydrophobe;
		out << " " << info.xs_donor;
		out << " " << info.xs_acceptor;
		out << " " << info.ad_heteroatom;
		out << "\n";
	}
}


void setup_user_gd(grid_dims& gd, std::ifstream& user_in)
{
	std::string line;
	size_t pLines = 3;
	std::vector<std::string> temp;
	fl center_x = 0, center_y = 0, center_z = 0, size_x = 0, size_y = 0,
			size_z = 0;

	for (; pLines > 0; --pLines) //Eat first 3 lines
		std::getline(user_in, line);
	pLines = 3;

//Read in SPACING
	std::getline(user_in, line);
	boost::algorithm::split(temp, line, boost::algorithm::is_space());
	const fl granularity = ::atof(temp[1].c_str());
//Read in NELEMENTS
	std::getline(user_in, line);
	boost::algorithm::split(temp, line, boost::algorithm::is_space());
	size_x = (::atof(temp[1].c_str()) + 1) * granularity; // + 1 here?
	size_y = (::atof(temp[2].c_str()) + 1) * granularity;
	size_z = (::atof(temp[3].c_str()) + 1) * granularity;
//Read in CENTER
	std::getline(user_in, line);
	boost::algorithm::split(temp, line, boost::algorithm::is_space());
	center_x = ::atof(temp[1].c_str()) + 0.5 * granularity;
	center_y = ::atof(temp[2].c_str()) + 0.5 * granularity;
	center_z = ::atof(temp[3].c_str()) + 0.5 * granularity;

	vec span(size_x, size_y, size_z);
	vec center(center_x, center_y, center_z);
	VINA_FOR_IN(i, gd)
	{
		gd[i].n = sz(std::ceil(span[i] / granularity));
		fl real_span = granularity * gd[i].n;
		gd[i].begin = center[i] - real_span / 2;
		gd[i].end = gd[i].begin + real_span;
	}

}

//enum options and their parsers
enum ApproxType
{
	LinearApprox, SplineApprox, Exact, GPU
};

std::istream& operator>>(std::istream& in, ApproxType& type)
{
	using namespace boost::program_options;

	std::string token;
	in >> token;
	if (token == "spline")
		type = SplineApprox;
	else if (token == "linear")
		type = LinearApprox;
	else if (token == "exact")
		type = Exact;
	else if (token == "gpu")
	  type = GPU;
	else
		throw validation_error(validation_error::invalid_option_value);
	return in;
}

//set the default device to device and exit with error if there are any problems
static void initializeCUDA(int device)
{
	cudaSetDevice(device);
	cudaError_t error;
	cudaDeviceProp deviceProp;
	error = cudaGetDevice(&device);

	if (error != cudaSuccess)
	{
		std::cerr << "cudaGetDevice returned error code " << error << "\n";
		exit(-1);
	}

	error = cudaGetDeviceProperties(&deviceProp, device);

	if (deviceProp.computeMode == cudaComputeModeProhibited)
	{
		std::cerr << "Error: device is running in <Compute Mode Prohibited>, no threads can use ::cudaSetDevice().\n";
		exit(-1);
	}

	if (error != cudaSuccess)
	{
		std::cerr << "cudaGetDeviceProperties returned error code " << error << "\n";
		exit(-1);
	}
}


int main(int argc, char* argv[])
{
	using namespace boost::program_options;
	const std::string version_string =
			"Smina " __DATE__ ".  Based on AutoDock Vina 1.1.2.";
	const std::string error_message =
        "\n\n Please report this error at http://smina.sf.net\n"
        "Please remember to include the following in your problem report:\n"
        "* the EXACT error message,\n"
        "* your version of the program,\n"
        "* the type of computer system you are running it on,\n"
        "* all command line options,\n"
        "* configuration file (if used),\n"
        "* ligand file as PDBQT,\n"
        "* receptor file as PDBQT,\n"
        "* flexible side chains file as PDBQT (if used),\n"
        "* output file as PDBQT (if any),\n"
        "* input (if possible),\n"
        "* random seed the program used (this is printed when the program starts).\n"
        "\n Thank you!\n";

	const std::string cite_message =
"              _             \n" \
"             (_)            \n" \
"   __ _ _ __  _ _ __   __ _ \n" \
"  / _` | '_ \\| | '_ \\ / _` |\n" \
" | (_| | | | | | | | | (_| |\n" \
"  \\__, |_| |_|_|_| |_|\\__,_|\n" \
"   __/ |                    \n"  \
"  |___/                     \n"  \
					"\ngnina is based off of smina and AutoDock Vina.\nPlease cite appropriately.\n\n" \
"*** IMPORTANT: gnina is not yet intended for production use. Use smina. ***\n";

	try
	{
		std::string rigid_name, flex_name, config_name, log_name, atom_name;
		std::vector<std::string> ligand_names;
		std::string out_name;
		std::string outf_name;
		std::string ligand_names_file;
		std::string atomconstants_file;
		std::string custom_file_name;
		std::string usergrid_file_name;
		std::string flex_res;
		double flex_dist = -1.0;
		fl center_x = 0, center_y = 0, center_z = 0, size_x = 0, size_y = 0,
				size_z = 0;
		fl autobox_add = 4;
		fl out_min_rmsd = 1;
		std::string autobox_ligand;
		std::string flexdist_ligand;
		std::string builtin_scoring;
		int device = 0;

		// -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 0.05846
		fl weight_gauss1 = -0.035579;
		fl weight_gauss2 = -0.005156;
		fl weight_repulsion = 0.840245;
		fl weight_hydrophobic = -0.035069;
		fl weight_hydrogen = -0.587439;
		fl weight_rot = 0.05846;
		fl user_grid_lambda;
		bool help = false, help_hidden = false, version = false;
		bool quiet = false;
		bool accurate_line = false;
		bool flex_hydrogens = false;
		bool gpu_on = false;
		bool print_terms = false;
		bool print_atom_types = false;
		bool add_hydrogens = true;
		bool no_lig = false;

		user_settings settings;
		minimization_params minparms;
		ApproxType approx = LinearApprox;
		fl approx_factor = 32;

		positional_options_description positional; // remains empty

		options_description inputs("Input");
		inputs.add_options()
		("receptor,r", value<std::string>(&rigid_name),
				"rigid part of the receptor")
		("flex", value<std::string>(&flex_name),
				"flexible side chains, if any (PDBQT)")
		("ligand,l", value<std::vector<std::string> >(&ligand_names),
				"ligand(s)")
		("flexres", value<std::string>(&flex_res),
				"flexible side chains specified by comma separated list of chain:resid")
		("flexdist_ligand", value<std::string>(&flexdist_ligand),
						"Ligand to use for flexdist")
		("flexdist", value<double>(&flex_dist),
				"set all side chains within specified distance to flexdist_ligand to flexible");

		//options_description search_area("Search area (required, except with --score_only)");
		options_description search_area("Search space (required)");
		search_area.add_options()
		("center_x", value<fl>(&center_x), "X coordinate of the center")
		("center_y", value<fl>(&center_y), "Y coordinate of the center")
		("center_z", value<fl>(&center_z), "Z coordinate of the center")
		("size_x", value<fl>(&size_x), "size in the X dimension (Angstroms)")
		("size_y", value<fl>(&size_y), "size in the Y dimension (Angstroms)")
		("size_z", value<fl>(&size_z), "size in the Z dimension (Angstroms)")
		("autobox_ligand", value<std::string>(&autobox_ligand),
				"Ligand to use for autobox")
		("autobox_add", value<fl>(&autobox_add),
				"Amount of buffer space to add to auto-generated box (default +4 on all six sides)")
		("no_lig", bool_switch(&no_lig)->default_value(false),
				"no ligand; for sampling/minimizing flexible residues");

		//options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
		options_description outputs("Output (optional)");
		outputs.add_options()
		("out,o", value<std::string>(&out_name),
				"output file name, format taken from file extension")
		("out_flex",value<std::string>(&outf_name),
				"output file for flexible receptor residues")
		("log", value<std::string>(&log_name), "optionally, write log file")
		("atom_terms", value<std::string>(&atom_name),
				"optionally write per-atom interaction term values")
		("atom_term_data", bool_switch(&settings.include_atom_info)->default_value(false),
				"embedded per-atom interaction terms in output sd data");

		options_description scoremin("Scoring and minimization options");
		scoremin.add_options()
		("scoring", value<std::string>(&builtin_scoring),"specify alternative builtin scoring function")
		("custom_scoring", value<std::string>(&custom_file_name),
				"custom scoring function file")
		("custom_atoms", value<std::string>(&atomconstants_file), "custom atom type parameters file")
		("score_only", bool_switch(&settings.score_only)->default_value(false), "score provided ligand pose")
		("local_only", bool_switch(&settings.local_only)->default_value(false),
				"local search only using autobox (you probably want to use --minimize)")
		("minimize", bool_switch(&settings.dominimize)->default_value(false),
				"energy minimization")
		("randomize_only", bool_switch(&settings.randomize_only),
				"generate random poses, attempting to avoid clashes")
		("minimize_iters",
				value<unsigned>(&minparms.maxiters)->default_value(0),
				"number iterations of steepest descent; default scales with rotors and usually isn't sufficient for convergence")
		("accurate_line", bool_switch(&accurate_line),
				"use accurate line search")
		("minimize_early_term", bool_switch(&minparms.early_term),
				"Stop minimization before convergence conditions are fully met.")
		("approximation", value<ApproxType>(&approx),
				"approximation (linear, spline, or exact) to use")
		("factor", value<fl>(&approx_factor),
				"approximation factor: higher results in a finer-grained approximation")
		("force_cap",value<fl>(&settings.forcecap),"max allowed force; lower values more gently minimize clashing structures")
		("user_grid", value<std::string>(&usergrid_file_name),
				"Autodock map file for user grid data based calculations")
		("user_grid_lambda", value<fl>(&user_grid_lambda)->default_value(-1.0),
								"Scales user_grid and functional scoring")
		("print_terms", bool_switch(&print_terms),
				"Print all available terms with default parameterizations")
				("print_atom_types", bool_switch(&print_atom_types), "Print all available atom types");

		options_description hidden("Hidden options for internal testing");
		hidden.add_options()
		("verbosity", value<int>(&settings.verbosity)->default_value(1),
				"Adjust the verbosity of the output, default: 1")
    ("flex_hydrogens", bool_switch(&flex_hydrogens),
        "Enable torsions effecting only hydrogens (e.g. OH groups). This is stupid but provides compatibility with Vina.");


		options_description misc("Misc (optional)");
		misc.add_options()
		("cpu", value<int>(&settings.cpu),
				"the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
		("seed", value<int>(&settings.seed), "explicit random seed")
		("exhaustiveness", value<int>(&settings.exhaustiveness)->default_value(8),
				"exhaustiveness of the global search (roughly proportional to time)")
		("num_modes", value<sz>(&settings.num_modes)->default_value(9),
				"maximum number of binding modes to generate")
		("energy_range", value<fl>(&settings.energy_range)->default_value(3.0),
				"maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
		("min_rmsd_filter", value<fl>(&out_min_rmsd)->default_value(1.0),
				"rmsd value used to filter final poses to remove redundancy")
		("quiet,q", bool_switch(&quiet), "Suppress output messages")
		("addH", value<bool>(&add_hydrogens),
				"automatically add hydrogens in ligands (on by default)")
		("device", value<int>(&device)->default_value(0), "GPU device to use")
		("gpu", bool_switch(&gpu_on), "Turn on GPU acceleration");

		options_description config("Configuration file (optional)");
		config.add_options()("config", value<std::string>(&config_name),
				"the above options can be put here");
		options_description info("Information (optional)");
		info.add_options()
		("help", bool_switch(&help), "display usage summary")
		("help_hidden", bool_switch(&help_hidden),
				"display usage summary with hidden options")
		("version", bool_switch(&version), "display program version");

		options_description desc, desc_simple;
		desc.add(inputs).add(search_area).add(outputs).add(scoremin).
				add(hidden).add(misc).add(config).add(info);
		desc_simple.add(inputs).add(search_area).add(scoremin).
				add(outputs).add(misc).add(config).add(info);

		variables_map vm;
		try
		{
			store(
					command_line_parser(argc, argv).options(desc)
							.style(
							command_line_style::default_style
									^ command_line_style::allow_guessing)
							.positional(positional).run(), vm);
			notify(vm);
		} catch (boost::program_options::error& e)
		{
			std::cerr << "Command line parse error: " << e.what() << '\n'
					<< "\nCorrect usage:\n" << desc_simple << '\n';
			return 1;
		}
		if (vm.count("config"))
		{
			try
			{
				ifile config_stream(config_name);
				store(parse_config_file(config_stream, desc), vm);
				notify(vm);
			} catch (boost::program_options::error& e)
			{
				std::cerr << "Configuration file parse error: " << e.what()
						<< '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
				return 1;
			}
		}
		if (help)
		{
			std::cout << desc_simple << '\n';
			return 0;
		}
		if (help_hidden)
		{
			std::cout << desc << '\n';
			return 0;
		}
		if (version)
		{
			std::cout << version_string << '\n';
			return 0;
		}

		if (!atomconstants_file.empty())
			setup_atomconstants_from_file(atomconstants_file);

		if (print_terms)
		{
			custom_terms t;
			t.print_available_terms(std::cout);
			return 0;
		}


		if(print_atom_types)
		{
			print_atom_info(std::cout);
			return 0;
		}

		initializeCUDA(device);

		set_fixed_rotable_hydrogens(!flex_hydrogens);

		if (settings.dominimize) //set default settings for minimization
		{
			if(!vm.count("force_cap"))
				settings.forcecap = 10; //nice and soft
			if (minparms.maxiters == 0)
				minparms.maxiters = 10000; //will presumably converge
			settings.local_only = true;
			minparms.type = minimization_params::BFGSAccurateLineSearch;

			if (!vm.count("approximation"))
				approx = SplineApprox;
			if (!vm.count("factor"))
				approx_factor = 10;
		}

		if (accurate_line)
		{
			minparms.type = minimization_params::BFGSAccurateLineSearch;
		}

		bool search_box_needed = !(settings.score_only || settings.local_only); // randomize_only and local_only still need the search space; dkoes - for local get box from ligand
		bool output_produced = !settings.score_only;
		bool receptor_needed = !settings.randomize_only;

		if (receptor_needed)
		{
			if (vm.count("receptor") <= 0)
			{
				std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n"
						<< desc_simple << '\n';
				return 1;
			}
		}

		if (ligand_names.size() == 0)
		{
			if (!no_lig)
			{
				std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n"
						<< desc_simple << '\n';
				return 1;
			}
			else //put in "fake" ligand
			{
				ligand_names.push_back("");
			}
		}
		else if (no_lig) //ligand specified with no_lig
		{
			std::cerr << "Ligand specified with --no_lig.\n"
					<< "\nCorrect usage:\n"
					<< desc_simple << '\n';
			return 1;
		}

		if (settings.exhaustiveness < 1)
			throw usage_error("exhaustiveness must be 1 or greater");
		if (settings.num_modes < 1)
			throw usage_error("num_modes must be 1 or greater");

		boost::optional<std::string> flex_name_opt;
		if (vm.count("flex"))
			flex_name_opt = flex_name;

		if (vm.count("flex") && !vm.count("receptor"))
			throw usage_error(
					"Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

		tee log(quiet);
		if (vm.count("log") > 0)
			log.init(log_name);

		std::ofstream atomoutfile;
		if (vm.count("atom_terms") > 0)
			atomoutfile.open(atom_name.c_str());

		if (autobox_ligand.length() > 0)
		{
			setup_autobox(autobox_ligand, autobox_add,
					center_x, center_y, center_z,
					size_x, size_y, size_z);
		}

		if (search_box_needed && autobox_ligand.length() == 0)
		{
			options_occurrence oo = get_occurrence(vm, search_area);
			if (!oo.all)
			{
				check_occurrence(vm, search_area);
				std::cerr << "\nCorrect usage:\n" << desc_simple << std::endl;
				return 1;
			}
			if (size_x <= 0 || size_y <= 0 || size_z <= 0)
				throw usage_error("Search space dimensions should be positive");
		}

		if(flex_dist > 0 && flexdist_ligand.size() == 0)
		{
			throw usage_error("Must specify flexdist_ligand with flex_dist");
		}

		FlexInfo finfo(flex_res, flex_dist, flexdist_ligand, log);

		log << cite_message << '\n';

		grid_dims gd; // n's = 0 via default c'tor
		grid_dims user_gd;
		grid user_grid;

		flv weights;

		//dkoes, set the scoring function
		custom_terms t;
		if(user_grid_lambda != -1.0){
			t.set_scaling_factor(user_grid_lambda);
		}
		if (custom_file_name.size() > 0)
		{
			ifile custom_file(custom_file_name);
			t.add_terms_from_file(custom_file);
		}
		else if (builtin_scoring.size() > 0)
		{
			if(!builtin_scoring_functions.set(t, builtin_scoring))
			{
				std::stringstream ss;
				builtin_scoring_functions.print_functions(ss);
				throw usage_error("Invalid builtin scoring function: "+builtin_scoring+". Options are:\n"+ss.str());
			}
		}
		else
		{
			t.add("gauss(o=0,_w=0.5,_c=8)", -0.035579);
			t.add("gauss(o=3,_w=2,_c=8)", -0.005156);
			t.add("repulsion(o=0,_c=8)", 0.840245);
			t.add("hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069);
			t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439);
			t.add("num_tors_div", 5 * 0.05846 / 0.1 - 1);
		}
				
		log << std::setw(12) << std::left << "Weights" << " Terms\n" << t
				<< "\n";

		if (usergrid_file_name.size() > 0)
		{
			ifile user_in(usergrid_file_name);
			fl ug_scaling_factor = 1.0;
			if(user_grid_lambda != -1.0){
				ug_scaling_factor = 1 - user_grid_lambda;
			}
			setup_user_gd(user_gd, user_in);
			user_grid.init(user_gd, user_in, ug_scaling_factor); //initialize user grid
		}

		const fl granularity = 0.375;
		if (search_box_needed)
		{
			vec span(size_x, size_y, size_z);
			vec center(center_x, center_y, center_z);
			VINA_FOR_IN(i, gd)
			{
				gd[i].n = sz(std::ceil(span[i] / granularity));
				fl real_span = granularity * gd[i].n;
				gd[i].begin = center[i] - real_span / 2;
				gd[i].end = gd[i].begin + real_span;
			}
		}

		if (vm.count("cpu") == 0)
		{
			settings.cpu = boost::thread::hardware_concurrency();
			if (settings.verbosity > 1)
			{
				if (settings.cpu > 0)
					log << "Detected " << settings.cpu << " CPU"
							<< ((settings.cpu > 1) ? "s" : "") << '\n';
				else
					log << "Could not detect the number of CPUs, using 1\n";
			}
		}
		if (settings.cpu < 1)
			settings.cpu = 1;
		if (settings.verbosity > 1 && settings.exhaustiveness < settings.cpu)
			log
					<< "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";


		//dkoes - parse in receptor once
		MolGetter mols(rigid_name, flex_name, finfo, add_hydrogens, log);

		//dkoes, hoist precalculation outside of loop
		weighted_terms wt(&t, t.weights());

		boost::shared_ptr<precalculate> prec;

		if (gpu_on || approx == GPU)
		{ //don't get a choice
		  prec = boost::shared_ptr<precalculate>(new precalculate_gpu(wt, approx_factor));
		}
		else if (approx == SplineApprox)
			prec = boost::shared_ptr<precalculate>(
					new precalculate_splines(wt, approx_factor));
		else if (approx == LinearApprox)
			prec = boost::shared_ptr<precalculate>(
					new precalculate_linear(wt, approx_factor));
		else if (approx == Exact)
			prec = boost::shared_ptr<precalculate>(
					new precalculate_exact(wt));

		//setup single outfile
		using namespace OpenBabel;
		ozfile outfile;
		std::string outext;
		if (out_name.length() > 0)
		{
			outext = outfile.open(out_name);
		}

		ozfile outflex;
		std::string outfext;
		if(outf_name.length() > 0)
		{
			outfext = outflex.open(outf_name);
		}

		if(settings.score_only) //output header
		{
			std::vector<std::string> enabled_names = t.get_names(true);
			log << "## Name";
			VINA_FOR_IN(i, enabled_names)
			{
				log << " " << enabled_names[i];
			}
			for (unsigned i = 0, n = t.conf_independent_terms.size(); i < n; i++)
			{
				log << " " << t.conf_independent_terms[i].name;
			}
			log << "\n";
		}

	  boost::timer time;
		//loop over input ligands
		for (unsigned l = 0, nl = ligand_names.size(); l < nl; l++)
		{
			doing(settings.verbosity, "Reading input", log);
			const std::string& ligand_name = ligand_names[l];
			mols.setInputFile(ligand_name);

			//process input molecules one at a time
			unsigned i = 0;
            /* TODO */
            model *_m = new model;
            model &m = *_m;
			while (mols.readMoleculeIntoModel(m))
			{
				if (settings.local_only)
				{
					//dkoes - for convenience get box from model
					gd = m.movable_atoms_box(autobox_add, granularity);
				}

				boost::optional<model> ref;
				done(settings.verbosity, log);

				std::stringstream output;
				std::vector<result_info> results;

				main_procedure(m, *prec, ref, settings,
						false, // no_cache == false
						atomoutfile.is_open() || settings.include_atom_info, gpu_on,
						gd,minparms, wt,log, results, user_grid);

				if (outfile)
				{
					//write out molecular data
					for (unsigned j = 0, nr = results.size(); j < nr; j++)
					{
						results[j].write(outfile, outext, settings.include_atom_info, &wt, j+1);
					}
				}
				if(outflex)
				{
					//write out flexible residue data data
					for (unsigned j = 0, nr = results.size(); j < nr; j++)
					{
						results[j].writeFlex(outflex, outfext);
					}
				}
				if (atomoutfile)
				{
					for (unsigned j = 0, m = results.size(); j < m; j++)
					{
						results[j].writeAtomValues(atomoutfile, &wt);
					}
				}
				i++;
				if(no_lig) break; //only go through loop once

			}
		}

	  std::cout << "Loop time " << time.elapsed() << "\n";

	} catch (file_error& e)
	{
		std::cerr << "\n\nError: could not open \"" << e.name.string()
				<< "\" for " << (e.in ? "reading" : "writing") << ".\n";
		return 1;
	} catch (boost::filesystem::filesystem_error& e)
	{
		std::cerr << "\n\nFile system error: " << e.what() << '\n';
		return 1;
	} catch (usage_error& e)
	{
		std::cerr << "\n\nUsage error: " << e.what() << "\n";
		return 1;
	} catch (parse_error& e)
	{
		std::cerr << "\n\nParse error on line " << e.line << " in file \""
				<< e.file.string() << "\": " << e.reason << '\n';
		return 1;
	} catch (std::bad_alloc&)
	{
		std::cerr << "\n\nError: insufficient memory!\n";
		return 1;
	} catch (scoring_function_error e)
	{
		std::cerr << "\n\nError with scoring function specification.\n";
		std::cerr << e.msg << "[" << e.name << "]\n";
		return 1;
	}

// Errors that shouldn't happen:

	catch (std::exception& e)
	{
		std::cerr << "\n\nAn error occurred: " << e.what() << ". "
				<< error_message;
		return 1;
	} catch (internal_error& e)
	{
		std::cerr << "\n\nAn internal error occurred in " << e.file << "("
				<< e.line << "). " << error_message;
		return 1;
	} catch (...)
	{
		std::cerr << "\n\nAn unknown error occurred. " << error_message;
		return 1;
	}
}
