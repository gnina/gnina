#include <openbabel/mol.h>
#include <set>
#include "cnn_scorer.h"
#include "molgetter.h"

struct visualization_options
{
  bool frags_only;
  bool atoms_only;
  bool verbose;

  std::string outrec;
  std::string outlig;
};

class ColoredMol
{
    public:
    ColoredMol(std::string inLigName, std::string inRecName, std::string inModel, std::string inWeights, float inSize, std::string inOutRec, std::string inOutLig, const cnn_options &cnnopts, FlexInfo &finfo, tee &log, const vec &center, bool inNo_frag = false, bool inVerbose = false);
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
    bool no_frag, verbose;

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
