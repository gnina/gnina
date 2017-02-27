#include <iostream>
#include <iomanip>
#include <set>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/mol.h>
#include "cnn_visualization.hpp"
#include <boost/filesystem.hpp>
#include "cnn_scorer.h"
#include "molgetter.h"
#include "obmolopener.h"
#include "model.h"
#include "parse_pdbqt.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolOps.h>

using namespace OpenBabel;

cnn_visualization::cnn_visualization (const vis_options &visopts, const cnn_options &cnnopts, const vec &center )
{
    if(visopts.gpu > -1)
    {
        caffe::Caffe::SetDevice(visopts.gpu);
        caffe::Caffe::set_mode(caffe::Caffe::GPU);
    }

    this->visopts = visopts;
    this->cnnopts = cnnopts;
    this->center = &center;

    OBConversion conv;
    obmol_opener opener;

    try
    {
        opener.openForInput(conv, visopts.ligand_name);
        conv.Read(&lig_mol);
        opener.openForInput(conv, visopts.receptor_name);
        conv.Read(&rec_mol);
    }
    catch(file_error& e)
    {
        std::cout << "Could not open \"" << e.name.string() << "\" for reading\n";
        exit(1);
    }
}


void cnn_visualization::lrp()
{
    std::cout << "Starting LRP...\n";

    process_molecules();

    std::stringstream rec_stream(rec_string);
    model receptor = parse_receptor_pdbqt("", rec_stream);
    CNNScorer scorer(cnnopts, *center, receptor);

    std::stringstream lig_stream(lig_string);
    model ligand = parse_ligand_stream_pdbqt("", lig_stream);

    receptor.append(ligand);

    scorer.lrp(receptor);
}


void cnn_visualization::color()
{
    if(visopts.verbose)
    {
      print();
    }

    process_molecules();

    std::stringstream rec_stream(rec_string);
    unmodified_receptor = parse_receptor_pdbqt("", rec_stream);
    CNNScorer base_scorer(cnnopts, *center, unmodified_receptor);

    std::stringstream lig_stream(lig_string);
    unmodified_ligand = parse_ligand_stream_pdbqt("", lig_stream);
    model temp_rec = unmodified_receptor;

    temp_rec.append(unmodified_ligand);
    original_score = base_scorer.score(temp_rec);
    std::cout << "Original score: " << original_score << "\n\n";

    if(visopts.receptor_output.length() > 0)
    {
    remove_residues();
    }

    if(visopts.ligand_output.length() > 0)
    {
    remove_ligand_atoms();
    }
}

void cnn_visualization::print()
{
    std::cout << "ligand_name: " << visopts.ligand_name << '\n';
    std::cout << "receptor_name: " << visopts.receptor_name << '\n';
    std::cout << "cnn_model: " << cnnopts.cnn_model << '\n';
    std::cout << "cnn_weights: " << cnnopts.cnn_weights << '\n';
    std::cout << "box_size: " << visopts.box_size << '\n';
    std::cout << "receptor_output: " << visopts.receptor_output << '\n';
    std::cout << "ligand_output: " << visopts.ligand_output << '\n';
    std::cout << "frags_only: " << visopts.frags_only << '\n';
    std::cout << "atoms_only: " << visopts.atoms_only << '\n';
    std::cout << "verbose: " << visopts.verbose << "\n\n";
}

std::string cnn_visualization::modify_pdbqt(std::vector<int> atoms_to_remove, bool isRec)
{
    std::string mol_string;
    OBMol mol;

    std::sort(atoms_to_remove.begin(), atoms_to_remove.end()); //sort for efficient check when outputting pbdqt
    auto last = std::unique(atoms_to_remove.begin(), atoms_to_remove.end()); 
    atoms_to_remove.erase(last, atoms_to_remove.end()); //remove duplicates

    if(isRec)
    {
        mol_string = rec_string;
        mol = rec_mol;
    }
    else
    {
        mol_string = lig_string;
        mol = lig_mol;
    }

    int i = 0; //holds position in removal list

    if(visopts.verbose)
    {
    std::cout << "Removing: [";
    std::cout << atoms_to_remove[i];
    }

    std::stringstream ss;
    std::stringstream mol_stream(mol_string);

    ss << "ROOT\n"; //add necessary lines for gnina parsing

    bool list_ended = false;
    std::string line;
    while(std::getline(mol_stream, line))
    {
      if ((line.find("ATOM") < std::string::npos))
      {
        if(!list_ended)
        {
          std::string index_string = line.substr(7,5);
          int atom_index = std::stoi(index_string);
          if (atom_index != atoms_to_remove[i])
          {
            ss << line << '\n';
          }
          else
          {
            i++; //move to next item to check for skip
            if(i >= atoms_to_remove.size())
            {
              list_ended = true;
            }
            else if(visopts.verbose)
            {
            std::cout << ", " << atoms_to_remove[i];
            }
          }
        }
        else
        {
          ss << line << '\n';
        }
        
      }
    }
    ss << "ENDROOT\n";
    ss << "TORSDOF 0\n";

    if(visopts.verbose)
    {
    std::cout << "]\n";
    }

    return ss.str();
}

//add hydrogens with openbabel, generate PDBQT
//files for removal
void cnn_visualization::process_molecules()
{
    rec_mol.AddHydrogens();
    lig_mol.AddHydrogens(true, false, 7.4); //add only polar hydrogens

    OBConversion conv;

    conv.AddOption("r",OBConversion::OUTOPTIONS); //treat as rigid
    conv.AddOption("c",OBConversion::OUTOPTIONS); //combine rotatable portions of molecule
    conv.AddOption("p",OBConversion::OUTOPTIONS);
    conv.SetOutFormat("PDBQT"); //use pdbqt to make passing to parse_pdbqt possible

    //generate base ligand pdbqt string
    std::string temp_lig_string = conv.WriteString(&lig_mol);
    std::stringstream lig_stream;
    lig_stream << "ROOT\n";
    lig_stream << temp_lig_string;
    lig_stream << "ENDROOT\n" << "TORSDOF 0";
    lig_string = lig_stream.str();

    //generate base receptor pdbqt string
    std::string temp_rec_string = conv.WriteString(&rec_mol);
    std::stringstream rec_stream;
    rec_stream << "ROOT\n";
    rec_stream << temp_rec_string;
    rec_stream << "ENDROOT\n" << "TORSDOF 0";
    rec_string = rec_stream.str();

    OBMol lig_copy = lig_mol; //Center() will change atom coordinates, so use copy

    vector3 cen = lig_copy.Center(0);
    cenCoords[0] = cen.GetX();
    cenCoords[1] = cen.GetY();
    cenCoords[2] = cen.GetZ();
}

//scores provided receptor string against unmodified ligand
float cnn_visualization::score_modified_receptor(const std::string &modified_rec_string)
{
    std::stringstream rec_stream(modified_rec_string);
    std::stringstream lig_stream(lig_string);

    model m = parse_receptor_pdbqt("", rec_stream);

    CNNScorer cnn_scorer(cnnopts, *center, m);
    
    model l = parse_ligand_stream_pdbqt("", lig_stream);
    m.append(l);

    float score_val = cnn_scorer.score(m);
    if(visopts.verbose)
    {
    std::cout << "SCORE: " << score_val << '\n';
    }

    return score_val;
}

//scores provided ligand string against unmodified ligand
float cnn_visualization::score_modified_ligand(const std::string &mol_string)
{
    std::stringstream lig_stream(mol_string);
    std::stringstream rec_stream(rec_string);
    
    model temp = unmodified_receptor;
    static CNNScorer cnn_scorer; //share same scorer across all ligands

    static bool first = true;
    if(first)
    {
    cnn_scorer = CNNScorer(cnnopts, *center, temp);
    first = false;
    }

    model l = parse_ligand_stream_pdbqt("", lig_stream);
    temp.append(l);

    float score_val = cnn_scorer.score(temp);
    if(visopts.verbose)
    {
    std::cout << "SCORE: " << score_val << '\n';
    }

    return score_val;
}

//input: raw score differences
//will be transformed before writing
void cnn_visualization::write_scores(const std::vector<float> score_diffs, bool isRec)
{
    std::string file_name;
    std::string mol_string;

    if(isRec)
    {
        file_name = visopts.receptor_output;
        mol_string = rec_string;
    }
    else
    {
        file_name = visopts.ligand_output;
        mol_string = lig_string;
    }

    std::ofstream out_file; 
    out_file.open(file_name);

    out_file << "CNN MODEL: " << cnnopts.cnn_model << '\n';
    out_file << "CNN WEIGHTS: " << cnnopts.cnn_weights << '\n';

    std::stringstream mol_stream(mol_string);
    std::string line;
    std::string index_string;
    int atom_index;
    std::stringstream score_stream;
    std::string score_string;
    
    while(std::getline(mol_stream, line))
    {
        if ((line.find("ATOM") < std::string::npos) || 
            (line.find("HETATM") < std::string::npos))
        {
            score_stream.str(""); //clear stream for next score
            index_string = line.substr(6,5);
            atom_index = std::stoi(index_string);

            if ((score_diffs[atom_index] > 0.001) || (score_diffs[atom_index] < -0.001)) //ignore very small score differences
            {
                float transformed_score = transform_score_diff(score_diffs[atom_index]);
                score_stream << std::fixed << std::setprecision(5) << transformed_score;
                out_file << line.substr(0,61);
                score_string = score_stream.str();
                score_string.resize(5);
                out_file.width(5);
                out_file.fill('.');
                out_file << std::right << score_string;
                out_file << line.substr(66) << '\n';

            }
            else
            {
                out_file << line << '\n';
            }
        }
        else
        {
            out_file << line << '\n';
        }
    }
}


//returns false if not at least one atom within range of center of ligand
bool cnn_visualization::check_in_range(std::unordered_set<int> atoms)
{
    float x = cenCoords[0];
    float y = cenCoords[1];
    float z = cenCoords[2];

    float allowed_dist = visopts.box_size / 2;
    int num_atoms = rec_mol.NumAtoms();

    OBAtom* atom;

    for(auto i = atoms.begin(); i != atoms.end(); ++i)
    {

        if (*i >= num_atoms)
        {
            return false;
        }

        atom = rec_mol.GetAtom(*i);
        if(atom->GetX() < x + allowed_dist)
            if (atom->GetY() < y + allowed_dist)
                if (atom->GetZ() < z + allowed_dist)
                    if (atom->GetX() > x - allowed_dist)
                        if (atom->GetY() > y - allowed_dist)
                            if (atom->GetZ() > z - allowed_dist)
                                return true;
    }

    return false;
}

//transforms score diff to maximize digits in b-factor
//field, and to raise small values by square-rooting
float cnn_visualization::transform_score_diff(float diff_val)
{
    float temp = diff_val;

    if(temp < 0)
    {
        temp = 0 - std::sqrt(std::abs(temp));
    }
    else
    {
        temp = std::sqrt(temp);
    }

    temp = temp * 100;

    return temp;
}

//removes whole residues at a time, and scores the resulting receptor
void cnn_visualization::remove_residues()
{
    std::vector<float> score_diffs(rec_mol.NumAtoms() + 1, 0.00);
    std::string last_res = "";
    std::string curr_res;
    std::vector<int> atoms_to_remove;
    std::set<std::string> res_list;

    std::string mol_string = rec_string;
    std::stringstream mol_stream(mol_string);
    std::string line;

    //add all residues to res_list 
    while(std::getline(mol_stream, line))
    {
        if(line.find("ATOM") < std::string::npos)
        {
            curr_res = line.substr(23,4);
            if(line.substr(23,4) != last_res)
            {
                res_list.insert(curr_res);
                last_res = curr_res;
            }
        }
    }

    int res_count = res_list.size();
    int counter = 1;

    //iterate through residues
    for(auto i = res_list.begin(); i != res_list.end(); ++i)
    {
        if(!visopts.verbose)
        {
            std::cout << "Scoring residues: " << counter << '/' << res_count << '\r' << std::flush;
            counter++;
        }

        mol_stream.clear();
        mol_stream.str(mol_string);

        //add each atom in residue to removal list
        while(std::getline(mol_stream, line))
        {
            if ((line.find("ATOM") < std::string::npos))
            {
                if(line.substr(23,4).compare(*i) == 0)
                {
                    atoms_to_remove.push_back(std::stoi(line.substr(6,5)));
                }
            }
        }

        bool remove = true;

        if(!visopts.skip_bound_check)
        {
            //make set for check_in_range test
            std::unordered_set<int> remove_set;
            for(auto i: atoms_to_remove)
            {
                remove_set.insert(i);
            }

            remove = check_in_range(remove_set);
        }

        if(remove)
        {
            std::string modified_mol_string = modify_pdbqt(atoms_to_remove, true);

            float score_val = score_modified_receptor(modified_mol_string);

            for( auto f : atoms_to_remove)
            {
                if(f)
                {
                    float score_diff = original_score - score_val;
                    score_diffs[f] = score_diff;
                }
            }

            if(visopts.output_files)
            {
                output_modified_string(modified_mol_string, atoms_to_remove, true);
            }
        }

        atoms_to_remove.clear();
    }

    write_scores(score_diffs, true);
    std::cout << '\n';
}

//checks all input indices for hydrogen neighbors, and appends them
void cnn_visualization::add_adjacent_hydrogens(std::vector<int> &atoms_to_remove, bool isRec)
{
    OBMol mol;

    if(isRec)
    {
        mol = rec_mol;
        mol.AddHydrogens();
    }

    else
    {
        mol = lig_mol;
        mol.AddHydrogens();
    }

    int original_size = atoms_to_remove.size(); //atoms_to_remove will have atoms appended
    for(int i = 0; i < original_size; ++i)
    {
        OBAtom* atom = mol.GetAtom(atoms_to_remove[i]);
        OBAtom* neighbor;
        if (atoms_to_remove[i] > 0) //an index of 0 segfaults
        {
            for (OBAtomAtomIter neighbor(atom); neighbor; ++neighbor)
            {
                if(neighbor->GetAtomicNum() == 1)
                {
                    atoms_to_remove.push_back(neighbor->GetIdx());
                }
            }
        }
    }
}

//mostly for debug convenience
void cnn_visualization::print_vector(const std::vector<int> &atoms_to_remove)
{
    std::cout << "[" << atoms_to_remove[0];

    for(int i = 1; i < atoms_to_remove.size(); ++i)
    {
        std::cout << ", " << atoms_to_remove[i];
    }

    std::cout << "]";
}

//removes individual atoms, scores them, and returns the diffs
std::vector<float> cnn_visualization::remove_each_atom()
{
    OBMol lig_mol_h = lig_mol;
    lig_mol_h.AddHydrogens(); //just in case hydrogen numbers don't add up

    std::vector<float> score_diffs(lig_mol_h.NumAtoms(), 0);
    std::stringstream lig_stream(lig_string);
    std::string line;

    std::string index_string;
    int atom_index;
    std::vector<int> atoms_to_remove; 
    float score_val;

    int counter = 1;
    int num_atoms = lig_mol.NumAtoms();

    while (std::getline(lig_stream, line))
    {
        if (line.find("ATOM") < std::string::npos)
        {     
            if (!visopts.verbose)
            {
                std::cout << "Scoring individual atoms: " << counter << '/' << num_atoms << '\r' << std::flush;
                counter++;
            }

            std::string modified_mol_string;
            index_string = line.substr(6, 5);
            atom_index = std::stoi(index_string);

            if (lig_mol.GetAtom(atom_index)->GetAtomicNum() != 1) //don't remove hydrogens individually
            {
                atoms_to_remove.push_back(atom_index);
                add_adjacent_hydrogens(atoms_to_remove, false);

                modified_mol_string = modify_pdbqt(atoms_to_remove, false);
                score_val = score_modified_ligand(modified_mol_string);

                score_diffs[atom_index] = original_score - score_val;
            }

        if (visopts.output_files)
        {
            output_modified_string(modified_mol_string, atoms_to_remove, false);
        }
        
        atoms_to_remove.clear();

        }
    }

    if (visopts.verbose)
    {
        //print index:type for debugging
        for (auto i = lig_mol.BeginAtoms(); i != lig_mol.EndAtoms(); ++i)
        {
            std::cout << (*i)->GetIdx() << ": " << (*i)->GetType() << '\n';
        }
    }

    std::cout << '\n';
    return score_diffs;
}

//writes modified pdbqt strings to file, for debugging purposes
void cnn_visualization::output_modified_string(const std::string &modified_string, const std::vector<int> &atoms_removed,
        bool receptor)
{
    static bool first = true;
    static int rec_counter=0;
    static int lig_counter=0;

    if (first)
    {
        ofstream original_rec_out;
        original_rec_out.open("unmodified_receptor.pdbqt");
        original_rec_out << rec_string;
        original_rec_out.close();

        ofstream original_lig_out;
        original_lig_out.open("unmodified_ligand.pdbqt");
        original_lig_out << lig_string;
        original_lig_out.close();
      
        first = false;
    }

    ofstream file_out;
    std::stringstream filename;

    if (receptor)
    {
        filename << "mod_receptor_" << atoms_removed[0];
    }
    else
    {
        filename << "mod_ligand_" << atoms_removed[0];
    }

    filename << ".pdbqt";
    
    file_out.open(filename.str());

    file_out << "REMARK: ATOMS REMOVED [" << atoms_removed[0];

    for (int i = 1; i < atoms_removed.size(); ++i)
    {
        file_out << ", " << atoms_removed[i];
    }

    file_out << "]\n";
    file_out << modified_string;
    file_out.close();
}

//writes sum of score differences for each atom along with original score to
//file for analysis
void cnn_visualization::write_additivity(std::vector<float> single_score_diffs, std::vector<float> frag_score_diffs)
{
    float single_total = 0;
    float frag_total = 0;
    int num_atoms = lig_mol.NumAtoms();

    if(!visopts.frags_only)
    {
        for (int i = 1; i < single_score_diffs.size(); ++i)
        {
            if(i <= num_atoms)
            {
                if (lig_mol.GetAtom(i)->GetAtomicNum() != 1) //hydrogens will have score of 0
                {
                    single_total += single_score_diffs[i];
                }
            }
        }
    }

    if(!visopts.atoms_only)
    {
        for (int i = 1; i < frag_score_diffs.size(); ++i)
        {
            if(i <= num_atoms)
            {
                if(lig_mol.GetAtom(i)->GetAtomicNum() != 1) //hydrogens will have score of 0
                {
                    frag_total += frag_score_diffs[i];
                }
            }
        }
    }

    std::ofstream out_file;
    out_file.open(visopts.additivity, std::ios_base::app);

    if (visopts.verbose)
    {
            std::cout << "ORIGINAL SCORE: " << original_score << '\n';

            if(!visopts.frags_only)
            {
                std::cout << "SUM OF SINGLE REMOVALS: " << single_total << '\n';
            }
            if(!visopts.atoms_only)
            {
                std::cout << "SUM OF FRAGMENT REMOVALS: " << frag_total << '\n';
            }
    }

  out_file << visopts.ligand_name  << " " << original_score \
      << " " << single_total \
      << " " << frag_total \
      << "\n";

}

//wrapper for fragment and individual removals
void cnn_visualization::remove_ligand_atoms()
{
    std::vector <float> individual_score_diffs;
    std::vector <float> frag_score_diffs;
    if (visopts.atoms_only)
    {
        individual_score_diffs = remove_each_atom();
        write_scores(individual_score_diffs, false);
    }

    else if (visopts.frags_only)
    {
        frag_score_diffs = remove_fragments(6);
        write_scores(frag_score_diffs, false);
    }

    else
    {
        individual_score_diffs = remove_each_atom();
        frag_score_diffs = remove_fragments(6);

        std::vector<float> both_score_diffs (individual_score_diffs.size(), 0.00);

        //average individual and fragment diffs
        for(int i = 0; i != individual_score_diffs.size(); ++i)
        {
            float avg = (individual_score_diffs[i] + frag_score_diffs[i]) / 2;
            both_score_diffs[i] = avg;
        }

        write_scores(both_score_diffs, false);

    }

    if(visopts.additivity.length() > 0)
        {
            write_additivity(individual_score_diffs, frag_score_diffs);
        }


}

//returns openbabel index of atom with supplied x coordinate
//necessary because rdkit indexing seems inconsistent with openbabel
int cnn_visualization::get_openbabel_index(double x_coordinate)
{
    static bool first = true;
    static std::map<double, int> coord_to_index;

    //fill map on first run
    if (first)
    {
        for (int i = 1; i < lig_mol.NumAtoms(); ++i)
        {
            double * coord = lig_mol.GetAtom(i)->GetCoordinate();

            if (visopts.verbose)
            {
                std::cout << *coord << ":" << i << '\n';
            }

            coord_to_index[*coord] = i;
        }

        first = false;
    }

    return coord_to_index[x_coordinate];

}

//returns average score difference for each index across all fragments
std::vector<float> cnn_visualization::remove_fragments(int size)
{
    OBConversion conv;

    OBMol lig_mol_h = lig_mol;
    lig_mol_h.AddHydrogens(); //just in case hydrogen numbers don't add up

    std::vector<float> score_diffs(lig_mol_h.NumAtoms() + 1, 0.00);
    std::vector<int> score_counts(lig_mol_h.NumAtoms() + 1, 0); //used to calculate average of scores across all fragments
    
    //PDB has parsing issues with RDKit
    conv.SetOutFormat("MOL");
    std::stringstream MOL_stream(conv.WriteString(&lig_mol));

    unsigned int line = 0;

    RDKit::RWMol rdkit_mol(*(RDKit::MolDataStreamToMol(MOL_stream, line, false, true, false))); //removeHs = true
    RDKit::MolOps::removeHs(rdkit_mol, false, false, false); //hydrogens will be added by add_adjacent_hydrogens later

    if (visopts.verbose)
    {
        //print all bonds in rdkit_mol, to check against fragments
        for (auto i = rdkit_mol.beginBonds(); i != rdkit_mol.endBonds(); ++i)
        {
            int first = (*i)->getBeginAtomIdx() + 1;
            int second = (*i)->getEndAtomIdx() + 1;
            std::cout << first << "-" << second << '\n';
        }
    }

    std::vector<int> atoms_to_remove;

    //map of path length: list of paths
    RDKit::INT_PATH_LIST_MAP paths = RDKit::findAllSubgraphsOfLengthsMtoN(rdkit_mol, 1, size);

    RDKit::Conformer conf = rdkit_mol.getConformer();
 
    int path_count = 0;
    
    //count number of fragments for progress output
    for (auto path = paths.begin(); path != paths.end(); ++path) //iterate through path lengths
    {
        std::list<std::vector<int>> list_of_lists = std::get<1>(*path);
        for (auto bonds = list_of_lists.begin(); bonds != list_of_lists.end(); ++bonds) //iterate through paths of a given length
        {
            path_count++;
        }
    }
       
    int counter = 1; //stores current fragment number for progress output

    for (auto path = paths.begin(); path != paths.end(); ++path) //iterate through path lengths
    {
        std::list<std::vector<int>> list_of_lists = std::get<1>(*path);

        for (auto bonds = list_of_lists.begin(); bonds != list_of_lists.end(); ++bonds) //iterate through paths of a given length
        {
            std::vector<int> bond_list = *bonds;

            if (!visopts.verbose)
            {
                std::cout << "Scoring fragments: " << counter << '/' << path_count << '\r' << std::flush;
                counter++;
            }

            for (int i = 0; i < bond_list.size(); ++i) //iterate through bonds in path
            {
                RDKit::Bond bond = *(rdkit_mol.getBondWithIdx(bond_list[i]));
                int first = bond.getBeginAtomIdx();
                int second = bond.getEndAtomIdx();
                double first_coord = conf.getAtomPos(first).x;
                double second_coord = conf.getAtomPos(second).x;
                atoms_to_remove.push_back(get_openbabel_index(first_coord));
                atoms_to_remove.push_back(get_openbabel_index(second_coord));
            }

            int size_without_hydrogens = atoms_to_remove.size();
            add_adjacent_hydrogens(atoms_to_remove, false);

            std::string modified_ligand = modify_pdbqt(atoms_to_remove, false);
            float score = score_modified_ligand(modified_ligand);

            for(int i = 0; i < atoms_to_remove.size(); ++i)
            {
                int atom_index = atoms_to_remove[i];
                score_diffs[atoms_to_remove[i]] += (original_score - score) / size_without_hydrogens; //give each atom in removal equal portion of score difference
                score_counts[atoms_to_remove[i]] += 1;
            }

            atoms_to_remove.clear(); //clear for next group of atoms to be removed
        }
    }

    std::vector<float> avg_score_diffs(lig_mol.NumAtoms(), 0.00); 
    for (auto i = rdkit_mol.beginAtoms(); i != rdkit_mol.endAtoms(); ++i)
    {
        int r_index = (*i)->getIdx();
        int index = r_index + 1;
        if (score_counts[index] > 0)
        {  
            avg_score_diffs[index] = score_diffs[index] / score_counts[index];
            if (visopts.verbose)
            {
                std::cout << "Symbol: " << (*i)->getSymbol() << '\n';
                double x = (rdkit_mol.getConformer().getAtomPos(r_index)).x;
                std::cout << "X: " << x << '\n';
                std::cout << "RDKit Index: " << r_index << '\n';
                std::cout << "Corrected Index: " << index << '\n';
                std::cout << "Agg. Score Diff: " << score_diffs[index] << '\n';
                std::cout << "Score count: " << score_counts[index] << '\n';
                std::cout << "Avg. Score Diff: " << avg_score_diffs[index] << '\n';
                std::cout << "===============" << '\n';
            }
        }

        else
        {
            avg_score_diffs[index] = 0.00;
        }
    }

    std::cout << '\n';
    return avg_score_diffs;
    }

