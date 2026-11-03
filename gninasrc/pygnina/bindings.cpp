/* @THIS_IS_THE_SOURCE_FILE@
 * bindings.cpp
 *
 *  Python bindings for libmolgrid
 */

#include <boost/algorithm/string/trim.hpp>
#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <openbabel/obconversion.h>
#include <torch/torch.h>

#include <iostream>
#include <memory>
#include <sstream>
#include <streambuf>

#include "GninaConverter.h"
#include "cnn_torch_scorer.h"
#include "conf.h"
#include "custom_terms.h"
#include "dl_scorer.h"
#include "model.h"
#include "naive_non_cache.h"
#include "non_cache.h"
#include "parse_pdbqt.h"
#include "parsing.h"
#include "precalculate.h"
#include "python_streambuf.h"
#include "quasi_newton.h"
#include "result_info.h"
#include "user_opts.h"
#include "weighted_terms.h"

using namespace boost::python;
using namespace OpenBabel;
namespace bap = boost_adaptbx::python;

/** Python interface to gnina.  Maintains state (receptor file). */
class GNINA {
  model initm; // model with receptor only
  std::shared_ptr<DLScorer> cnn;
  minimization_params minparms;
  custom_terms t;
  std::shared_ptr<weighted_terms> wt;
  std::shared_ptr<precalculate> prec;
  std::shared_ptr<precalculate_exact> exact_prec;
  std::shared_ptr<naive_non_cache> nnc; // for scoring

  bool add_hydrogens = true;
  bool strip_hydrogens = false;

  GNINA(const GNINA &) = delete;
  GNINA &operator=(const GNINA &) = delete;

  void setLigand(std::istream &lig, const std::string &fmt, model &m) {
    if (initm.get_fixed_atoms().size() == 0) {
      throw usage_error("No receptor specified during scoring.");
    }
    m = initm;

    OBConversion conv;
    conv.SetOutFormat("PDBQT");

    OBFormat *format = conv.FormatFromExt(fmt);
    if (!format || !conv.SetInFormat(format)) {
      throw file_error(path(fmt), true);
    }

    OBMol mol;
    conv.Read(&mol, &lig); // will return after first success
    std::string name = mol.GetTitle();
    mol.StripSalts();
    m.set_name(name);

    if (mol.NumAtoms() == 0) {
      throw usage_error("Empty molecule");
    }

    parsing_struct p;
    context c;
    unsigned torsdof = GninaConverter::convertParsing(mol, p, c, add_hydrogens);
    non_rigid_parsed nr;
    postprocess_ligand(nr, p, c, torsdof);
    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());

    pdbqt_initializer tmp;
    tmp.initialize_from_nrp(nr, c, true);
    tmp.initialize(nr.mobility_matrix());
    if (strip_hydrogens)
      tmp.m.strip_hydrogens();

    m.append(tmp.m);
  }

  void init() {
    cnn_options cnnopts; // defaults for now

#ifdef USE_METAL
    if (torch::mps::is_available())
      cnn = std::make_shared<CNNTorchScorer<true>>(cnnopts);
    else
      cnn = std::make_shared<CNNTorchScorer<false>>(cnnopts);
#else
    if (torch::cuda::is_available())
      cnn = std::make_shared<CNNTorchScorer<true>>(cnnopts);
    else
      cnn = std::make_shared<CNNTorchScorer<false>>(cnnopts);
#endif

    minparms.maxiters = 10000;
    minparms.type = minimization_params::BFGSAccurateLineSearch;

    // vina scoring function
    t.add_vina();

    wt = std::make_shared<weighted_terms>(&t, t.weights());
    prec = std::make_shared<precalculate_splines>(*wt, 10.0);
    exact_prec = std::make_shared<precalculate_exact>(*wt);
    nnc = std::make_shared<naive_non_cache>(exact_prec.get());
  }

public:
  GNINA() { init(); }

  /** Initalize GNINA object with receptor contents */
  GNINA(const std::string &receptor, const std::string &format) {
    setReceptor(receptor, format);
    init();
  }

  GNINA(std::istream &receptor, const std::string &format) {
    setReceptor(receptor, format);
    init();
  }

  void setReceptor(const std::string &receptor, const std::string &format) {
    std::istringstream rec(receptor);
    setReceptor(rec, format);
  }

  void setReceptor(std::istream &rec, const std::string &format) {
    // molgetter code assumes file names so basically have to duplicate it here :-(
    OBConversion conv;
    conv.SetOutFormat("PDBQT");
    conv.AddOption("r", OBConversion::OUTOPTIONS);
    conv.AddOption("c", OBConversion::OUTOPTIONS);
    OBFormat *fmt = conv.FormatFromExt(format);
    if (!fmt || !conv.SetInFormat(fmt)) {
      throw usage_error("Incompatible format: " + format);
    }
    OBMol mol;
    if (!conv.Read(&mol, &rec))
      throw usage_error("Could not read receptor.");

    FOR_ATOMS_OF_MOL(a, mol) {
      OBResidue *residue = a->GetResidue();
      if (residue && a->GetFormalCharge() == 1) {
        std::string aname = residue->GetAtomID(&*a);
        boost::trim(aname);
        if (aname == "NH1") {
          a->SetFormalCharge(0);
          a->SetImplicitHCount(2);
        }
      }
    }

    mol.SetChainsPerceived(true);
    mol.AddHydrogens(true);
    FOR_ATOMS_OF_MOL(a, mol) { a->GetPartialCharge(); }

    std::string recstr = conv.WriteString(&mol);
    std::stringstream recstream(recstr);

    initm = parse_receptor_pdbqt("receptor", recstream);
  }

  /** Score provided ligand pose. */
  result_info score(const std::string &ligand, const std::string &format) {
    std::istringstream lig(ligand);
    return score(lig, format);
  }

  result_info score(std::istream &ligand, const std::string &format) {
    model m;
    float cnnscore = 0, cnnaffinity = 0, cnnvariance = 0, loss = 0;
    const vec authentic_v(1000, 1000, 1000);
    grid user_grid;

    setLigand(ligand, format, m);

    conf c = m.get_initial_conf(false);

    float intramolecular_energy = m.eval_intramolecular(*exact_prec, authentic_v, c);
    float e = m.eval_adjusted(*wt, *exact_prec, *nnc, authentic_v, c, intramolecular_energy, user_grid);

    cnn->set_center_from_model(m);
    cnnscore = cnn->score(m, false, cnnaffinity, loss, cnnvariance);

    return result_info(e, cnnscore, cnnaffinity, cnnvariance, -1, m);
  }

  result_info minimize(const std::string &ligand, const std::string &format, float scale = 0) {
    std::istringstream lig(ligand);
    return minimize(lig, format, scale);
  }

  result_info minimize(std::istream &ligand, const std::string &format, float scale = 0) { abort(); }
};

fl energy = 0;
fl cnnscore = -1;
fl cnnaffinity = 0;
fl cnnvariance = 0;
fl rmsd = -1;
std::string molstr;
std::string flexstr;
std::string atominfo;
std::string name;
bool sdfvalid = false;
;
bool flexsdfvalid = false;

BOOST_PYTHON_MODULE(pygnina) {
  Py_Initialize();

  class_<result_info, std::shared_ptr<result_info>>("result_info", "GNINA Results")
      .def("energy", &result_info::getEnergy)
      .def("cnnscore", &result_info::getCNNScore)
      .def("cnnaffinity", &result_info::getCNNAffinity)
      .def(
          "write", +[](result_info &self, object py_out, const std::string &ext) {
            bap::streambuf sb(py_out);
            std::ostream os(&sb);
            self.write(os, ext, false);
            os.flush();
          });

  class_<GNINA, std::shared_ptr<GNINA>, boost::noncopyable>("GNINA", "GNINA Python Interace")
      .def<void (GNINA::*)(const std::string &receptor, const std::string &format)>("set_receptor",
                                                                                    &GNINA::setReceptor)
      .def("set_receptor",
          +[](GNINA& self, object src, const std::string& ext) {
              if (PyUnicode_Check(src.ptr()) || PyBytes_Check(src.ptr())) {
                  self.setReceptor(extract<std::string>(src), ext);
              } else {
                  bap::streambuf sb(src);
                  std::istream is(&sb);
                  self.setReceptor(is, ext);
              }
          },
          (arg("source"), arg("format")))
      .def("score",
          +[](GNINA& self, object src, const std::string& ext) -> result_info {
              if (PyUnicode_Check(src.ptr()) || PyBytes_Check(src.ptr())) {
                  return self.score(extract<std::string>(src), ext);
              } else {
                  bap::streambuf sb(src);
                  std::istream is(&sb);
                  return self.score(is, ext);
              }
          },
          (arg("source"), arg("format")));          
}
