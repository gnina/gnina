#include <openbabel/mol.h>
#include <set>
#include "cnn_scorer.h"
#include "molgetter.h"

struct vis_options
{
  std::string ligName;
  std::string recName;
  std::string outRec;
  std::string outLig;

  bool frags_only;
  bool atoms_only;
  bool verbose;

  float size;

  vis_options(): size(23.5), frags_only(false), atoms_only(false), verbose(false) {}
};

class ColoredMol
{
    public:
    ColoredMol(const vis_options &visopts, const cnn_options &cnnopts, FlexInfo &finfo, tee &log, const vec &center);
    void color();
    void print();

    private:
    std::string ligName, recName, hRec, hLig, recPDB, cnnmodel, weights, outRec, outLig; 
    OpenBabel::OBMol ligMol, recMol, hLigMol, hRecMol;
    float size;
    float cenCoords [3];
    cnn_options cnnopts;
    FlexInfo* finfo;
    tee* log;
    const vec* center;
    bool frags_only, atoms_only,  verbose;

    void addHydrogens();
    float removeAndScore(std::vector<bool> removeList, bool isRec);
    void ligCenter();
    float score();
    void writeScores(std::vector<float> scoreList, bool isRec);
    bool inRange(std::set<int> atomList);
    std::vector<float> transform(std::vector<float> inList);
    void removeResidues();
    void removeEachAtom();
};
