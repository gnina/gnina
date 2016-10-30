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

using namespace OpenBabel;

cnn_visualization::cnn_visualization (const vis_options &visopts, const cnn_options &cnnopts, FlexInfo &finfo, tee &log, const vec &center )
    {
      this->visopts = visopts;
      this->cnnopts = cnnopts;
      this->finfo = &finfo;
      this->log = &log;
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

void cnn_visualization::color()
{
    if(visopts.verbose)
    {
      print();
    }

    process_molecules();
    ligCenter();

    std::stringstream rec_stream(rec_string);
    unmodified_receptor = parse_receptor_pdbqt("", rec_stream);
    CNNScorer base_scorer(cnnopts, *center, unmodified_receptor);

    std::stringstream lig_stream(lig_string);
    unmodified_ligand = parse_ligand_stream_pdbqt("", lig_stream);

    unmodified_receptor.append(unmodified_ligand);
    original_score = base_scorer.score(unmodified_receptor, true);
    std::cout << "Original Score: " << original_score << "\n\n";

    
    if(visopts.receptor_output.length() > 0)
    {
    std::cout << "Scoring residue removals:\n\n";
    remove_residues();
    }

    if(visopts.ligand_output.length() > 0)
    {
    std::cout << "Scoring individual atom removals:\n\n";
    remove_each_atom();
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

    std::cout << "Removing: [";

    std::stringstream ss;
    std::stringstream mol_stream(mol_string);

    ss << "ROOT\n"; //add necessary lines for gnina parsing


    int* i = atoms_to_remove.data();
    std::cout << *i;
    bool list_ended = false;
    std::string line;
    while(std::getline(mol_stream, line))
    {
      if ((line.find("ATOM") < std::string::npos))
      {
        if(! list_ended )
        {
          std::string index_string = line.substr(7,5);
          int atom_index = std::stoi(index_string);
          if (atom_index != *i)
          {
            ss << line << '\n';
          }
          else
          {
            ++i; //move to next item to check for skip
            if(i == &atoms_to_remove.back())
            {
              list_ended = true;
            }
            else
            {
            std::cout << ", " << *i;
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

    std::cout << "]\n";

    return ss.str();

    //ouput modified molecules for debugging
    /*
    static int counter = 0;
    std::stringstream file_name;
    file_name << "lig" << counter << ".pdbqt";
    counter++;
    std::ofstream file_out;
    file_out.open(file_name.str());
    file_out << ss.str();
    file_out.close();
    return scoreVal;
    */
}

//add hydrogens with openbabel, store PDB files for output, generate PDBQT
//files for removal
void cnn_visualization::process_molecules()
{
    rec_mol.AddHydrogens();
    lig_mol.AddHydrogens();

    OBConversion conv;
    conv.SetOutFormat("PDB");

    rec_pdb_string = conv.WriteString(&rec_mol); //store pdb's for score output
    lig_pdb_string = conv.WriteString(&lig_mol); 

    conv.SetOutFormat("PDBQT"); //use pdbqt to make passing to parse_pdbqt possible
    conv.AddOption("r",OBConversion::OUTOPTIONS);
    conv.AddOption("c",OBConversion::OUTOPTIONS);

    std::string temp_lig_string = conv.WriteString(&lig_mol);
    std::stringstream lig_stream;
    lig_stream << "ROOT\n";
    lig_stream << temp_lig_string;
    lig_stream << "ENDROOT\n" << "TORSDOF 0";
    lig_string = lig_stream.str();

    std::string temp_rec_string = conv.WriteString(&rec_mol);
    std::stringstream rec_stream;
    rec_stream << "ROOT\n";
    rec_stream << temp_rec_string;
    rec_stream << "ENDROOT\n" << "TORSDOF 0";
    rec_string = rec_stream.str();


    //conv.SetInFormat("PDBQT");
    //conv.ReadString(&hRecMol, hRec);
    //conv.ReadString(&hLigMol, hLig);
}

float cnn_visualization::score_modified_receptor(const std::string &modified_rec_string)
{
    std::stringstream rec_stream(modified_rec_string);
    std::stringstream lig_stream(lig_string);

    model m = parse_receptor_pdbqt("", rec_stream);

    CNNScorer cnn_scorer(cnnopts, *center, m);
    
    model l = parse_ligand_stream_pdbqt("", lig_stream);
    m.append(l);

    float score_val = cnn_scorer.score(m, true);
    std::cout << "SCORE: " << score_val << '\n';

    return score_val;
    //return 1;
}

float cnn_visualization::score_modified_ligand(const std::string &mol_string)
{
    std::stringstream lig_stream(mol_string);
    std::stringstream rec_stream(rec_string);
    
    model m = parse_receptor_pdbqt("", rec_stream);
    CNNScorer cnn_scorer(cnnopts, *center, m);

    model l = parse_ligand_stream_pdbqt("", lig_stream);
    unmodified_receptor.append(l);

    m.append(l);

    //float score_val = base_scorer.score(unmodified_receptor);
    float score_val = cnn_scorer.score(m, true);
    std::cout << "SCORE: " << score_val << '\n';

    return score_val;
    //return 1;
}


void cnn_visualization::write_scores(std::vector<float> scores, bool isRec)
{
    std::string file_name;
    std::string mol_string;
    if(isRec)
    {
        file_name = visopts.receptor_output;
        mol_string = rec_pdb_string;
    }
    else
    {
        file_name = visopts.ligand_output;
        mol_string = lig_pdb_string;
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

            if ((scores[atom_index] > 0.001) || (scores[atom_index] < -0.001)) //ignore very small scores
            {
                score_stream << std::fixed << std::setprecision(5) << scores[atom_index];
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


bool cnn_visualization::check_in_range(std::unordered_set<int> atoms)
{
    float x = cenCoords[0];
    float y = cenCoords[1];
    float z = cenCoords[2];

    float allowed_dist = visopts.box_size / 2;
    int num_atoms = rec_mol.NumAtoms();

    OBAtom* atom;
    for( auto i = atoms.begin(); i != atoms.end(); ++i)
    {
        if (*i >= num_atoms)
        {
            std::cout << "in range out of atoms";
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

void cnn_visualization::ligCenter()
{
    vector3 cen = lig_mol.Center(0);
    cenCoords[0] = cen.GetX();
    cenCoords[1] = cen.GetY();
    cenCoords[2] = cen.GetZ();
}

//calculates score difference and transforms to maximize digits in b-factor
//field
float cnn_visualization::transform_score(float score_val)
{
    float temp;
    temp = original_score - score_val;

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
void cnn_visualization::remove_residues()
{
    std::vector<float> scores(rec_mol.NumAtoms() + 1, 0.00);
    std::string last_res = "";
    std::string curr_res;
    std::vector<int> atoms_to_remove;
    std::set<std::string> res_list;

    std::string mol_string = rec_string;
    std::stringstream mol_stream(mol_string);
    std::string line;
    while(std::getline(mol_stream, line))
    {
      if ((line.find("ATOM") < std::string::npos))
        {
            curr_res = line.substr(23,4);
            if(line.substr(23,4) != last_res)
            {
                res_list.insert(curr_res);
                last_res = curr_res;
            }
        }
    }

    for( auto i = res_list.begin(); i != res_list.end(); ++i)
    {
        mol_stream.clear();
        mol_stream.str(mol_string);
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

        std::unordered_set<int> remove_set; //make set for check_in_range test
        //for (int i = 0; i < atoms_to_remove.size(); ++i)
        for(auto i: atoms_to_remove)
        {
          remove_set.insert(i);
        }
        if (check_in_range(remove_set))
        {
          std::string modified_mol_string = modify_pdbqt(atoms_to_remove, true);

          float score_val = score_modified_receptor(modified_mol_string);

          for ( auto f : atoms_to_remove)
          {
              if(f)
              {
              scores[f] = transform_score(score_val);
              }
          }
        }
        

        atoms_to_remove.clear();
    }

    write_scores(scores, true);
}

void cnn_visualization::remove_each_atom()
{
    std::vector<float> scores(lig_mol.NumAtoms());
    std::stringstream lig_stream(lig_string);
    std::string line;

    std::string index_string;
    int atom_index;
    std::vector<int> atoms_to_remove; 
    float score_val;

    while(std::getline(lig_stream, line))
    {
        if ((line.find("ATOM") < std::string::npos))
        {
            index_string = line.substr(6, 5);
            atom_index = std::stoi(index_string);
            if (lig_mol.GetAtom(atom_index)->GetAtomicNum() != 1) //dont remove hydrogens individually
            {
                atoms_to_remove.push_back(atom_index);
                OBAtom* atom;
                atom = lig_mol.GetAtom(atom_index);
                FOR_NBORS_OF_ATOM(neighbor, atom) //add any adjacent hydrogens
                {
                    if(neighbor->GetAtomicNum() == 1)
                    {
                        atoms_to_remove.push_back(neighbor->GetIdx());
                    }
                }

                
                std::string modified_mol_string = modify_pdbqt(atoms_to_remove, false);
                score_val = score_modified_ligand(modified_mol_string);

                scores[atom_index] = transform_score(score_val);
            }
        }
        atoms_to_remove.clear();

    }

    write_scores(scores, false);
}



