/*
 * molgridder.h
 *
 *  Created on: Apr 23, 2019
 *      Author: dkoes
 */

#ifndef MOLGRIDDER_H_
#define MOLGRIDDER_H_

#include <libmolgrid/grid_io.h>
#include <libmolgrid/atom_typer.h>
#include <libmolgrid/managed_grid.h>
#include <libmolgrid/grid_maker.h>
#include <vector>
#include "molgetter.h"
#include "gridoptions.h"

/** Store molecular model and apply gridermaker as needed to write out grids.
 *
 */
class MolGridder {
    double dimension = 0;
    double resolution = 0;
    bool random_rotate = false;
    float random_translate = 0.0;
    unsigned N = 0; //number of points on each side
    std::shared_ptr<libmolgrid::AtomTyper> rectyper;
    std::shared_ptr<libmolgrid::AtomTyper> ligtyper;

    std::vector<libmolgrid::MGrid3f> usergrids;

    //next read in receptor
    MolGetter mols; //use gnina routines for reading molecule
    libmolgrid::MGrid4f grid;
    libmolgrid::GridMaker gmaker;
    libmolgrid::Transform current_transform;

    gfloat3 center { 0, 0, 0 };
    bool center_set = false;
    bool gpu = false; //use gpu

    libmolgrid::Example ex; //coordinate/type data
    //set receptor from model into example
    void setReceptor(const model& m);
    //set ligand into ex from model (overwrite current)
    void setLigand(const model& m);
    //set grid from example
    void setGrid(bool use_gpu);

    //sets grid on cpu and compares to current
    void cpuSetGridCheck();

  public:
    MolGridder(const gridoptions& opt);

    //fix a center (instead of using ligand)
    void set_center(gfloat3 c);

    //return true if using fixed center
    bool has_set_center() const { return center_set; }

    //set grid parameters from an example dx
    void set_from_example(const std::string& examplefile);

    //read in specified grid files and set usergrids and grid parameters appropriately
    void set_usergrids(const std::vector<std::string>& userfiles);

    //output map for each grid
    void outputMAP(const std::string& base);

    //output an dx map for each grid
    void outputDX(const std::string& base);

    //output binary form of raw data in 3D multi-channel form
    void outputBIN(const std::string& base, bool outputrec, bool outputlig);

    //read a molecule (return false if unsuccessful)
    //set the ligand grid appropriately
    bool readMolecule(bool timeit);
};


#endif /* MOLGRIDDER_H_ */
