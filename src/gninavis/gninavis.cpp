#include <iostream>
#include "../lib/visualize.hpp"
#include "../lib/cnn_scorer.h"
#include "../lib/tee.h"
#include "../lib/flexinfo.h"
#include "../lib/box.h"
#include <vector>

int main()
{
    //test input
    int x [5] = {1, 2, 3, 4, 5};
    std::string ligName = "3gvu_lig.pdb";
    std::string recName = "3gvu_rec.pdb";
    std::string model = "/home/josh/models/refmodel.model";
    std::string weights = "/home/josh/models/refmodel.0_iter_10000.caffemodel";
    float size = 23.5;
    std::string outRec = "colored_3gvu_rec.pdb";
    std::string outLig = "colored_3gvu_lig.pdb";

    //necessary cnn_options
    cnn_options cnnopts;
    cnnopts.cnn_model = model;
    cnnopts.cnn_weights = weights;
    cnnopts.cnn_rotations = 24;
    cnnopts.cnn_scoring = true;

    //to pass to ColoredMol
    vis_options visopts;
    visopts.ligName = ligName;
    visopts.recName = recName;
    visopts.outRec = outRec;
    visopts.outLig = outLig;
    
    bool quiet = false; //maybe make option later
    tee log(quiet); //have to pass for cnn scoring

    //placeholders for FlexInfo to instantiate 
    std::string flex_res = "";
    float flex_dist = -1.0;
    std::string flexdist_ligand = "";
    FlexInfo finfo(flex_res, flex_dist, flexdist_ligand, log);
    
    //placeholders for center to instantiate
    float center_x = 0, center_y = 0, center_z = 0;
    vec center(center_x,center_y,center_z);

    ColoredMol cMol = ColoredMol(visopts, cnnopts, finfo, log, center);
    cMol.print();
    cMol.color();

    return 0;
}
