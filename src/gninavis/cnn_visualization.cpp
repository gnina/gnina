#include <iostream>
#include <iomanip>
#include <set>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/mol.h>
#include "visualize.hpp"
#include <boost/filesystem.hpp>
#include "cnn_scorer.h"
#include "molgetter.h"
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
    }

void cnn_visualization::color()
{
    OBConversion conv;

    std::string ext = boost::filesystem::extension(visopts.ligand_name);
    if(ext.compare(".pdb") == 0)
    {
        conv.SetInFormat("PDB");
    }
    else if(ext.compare(".pdbqt") == 0)
    {
        conv.SetInFormat("PDBQT");
    }
    else
    {
        std::cout << "File extension not supported: " << visopts.ligand_name << '\n';
        std::cout << "Please use .pdb or .pdbqt for ligand\n";
        exit(0);
    }

    conv.ReadFile(&ligMol, visopts.ligand_name);

    ext = boost::filesystem::extension(visopts.receptor_name);
    if(ext.compare(".pdb") == 0)
    {
      conv.SetInFormat("PDB");
    }
    else
    {
        std::cout << "File extension not supported: " << visopts.receptor_name << '\n';
        std::cout << "Please use .pdb for receptor\n";
        exit(0);
    }

    conv.ReadFile(&recMol, visopts.receptor_name);

    process_molecules();
    ligCenter();

    
    remove_residues();
    remove_each_atom();
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
    std::cout << "verbose: " << visopts.verbose << '\n';
}

float cnn_visualization::remove_and_score(std::vector<bool> removeList, bool isRec)
{
    std::string molString;
    OBMol mol;
    if(isRec)
    {
        molString = hRec;
        mol = hRecMol;
    }
    else
    {
        molString = hLig;
        mol = hLigMol;
    }

    std::cout << "Removing: [";
    if(!(isRec)) //if ligand
    {
        OBAtom* atom;
        for(int i = 0;i < removeList.size(); ++i)
        {
            if (removeList[i]) //index is in removeList
            {
            std::cout << i << "| ";
            atom = mol.GetAtom(i);
            FOR_NBORS_OF_ATOM(neighbor, atom)
            {
                if(neighbor->GetAtomicNum() == 1)
                {
                    //std::cout << "adding: " << neighbor->GetIdx() << '\n';
                    removeList[neighbor->GetIdx()] = true;
                }
            }
            }
        }
    }
    else //if receptor
    {
        //make set for check_in_range test
        std::set<int> removeSet;
        for (int i = 0; i < removeList.size(); ++i)
        {
            if (removeList[i])
            {
                removeSet.insert(i);
                std::cout << i << "| ";
            }
        }

       if (!(check_in_range(removeSet)))
        {
            return 0.00;
        }
    }

    std::cout << "]\n";
    std::stringstream ss;
    std::stringstream molStream(molString);

    ss << "ROOT\n"; //add necessary lines for gnina parsing

    std::string line;
    while(std::getline(molStream, line))
    {
        if((line.find("HETATM") < std::string::npos) ||
            (line.find("ATOM") < std::string::npos))
        {
            std::string firstNumString = line.substr(7,5);
            int firstNum = std::stoi(firstNumString);

            if (!(removeList[firstNum])) //don't write line if in list
            {
                ss << line << '\n';
            }


        }
    }
    ss << "ENDROOT\n";
    ss << "TORSDOF 0\n";

    float scoreVal = score(ss.str(), false);
    static int counter = 0;
    std::stringstream fileName;
    fileName << "lig" << counter << ".pdbqt";
    counter++;
    std::ofstream fileOut;
    fileOut.open(fileName.str());
    fileOut << ss.str();
    fileOut.close();
    return scoreVal;
}

//add hydrogens with openbabel, store PDB files for output, generate PDBQT
//files for removal
void cnn_visualization::process_molecules()
{
    recMol.AddHydrogens();
    ligMol.AddHydrogens();

    OBConversion conv;
    conv.SetInFormat("PDB");

    conv.SetOutFormat("PDB");
    recPDB = conv.WriteString(&recMol); //store pdb's for score output
    ligPDB = conv.WriteString(&ligMol); 

    conv.SetOutFormat("PDBQT"); //use pdbqt to make passing to parse_pdbqt possible
    conv.AddOption("r",OBConversion::OUTOPTIONS);
    conv.AddOption("c",OBConversion::OUTOPTIONS);
    hLig = conv.WriteString(&ligMol);
    hRec = conv.WriteString(&recMol);

    conv.SetInFormat("PDBQT");
    conv.ReadString(&hRecMol, hRec);
    conv.ReadString(&hLigMol, hLig);
}

float cnn_visualization::score(const std::string &molString, bool isRec)
{
  if( !isRec )
  {
        std::stringstream ligStream(molString);
        std::stringstream recStream(hRec);

        model m = parse_receptor_pdbqt("", recStream);

        CNNScorer cnn_scorer(cnnopts, *center, m);
        
        model l = parse_ligand_stream_pdbqt(molString, ligStream);
        m.append(l);

        float theScore = cnn_scorer.score(m);
        std::cout << "SCORE: " << theScore << '\n';

        return theScore;
  }
  else
  {
    std::cout << "receptor not implemented, exiting\n";
    exit(0);
  }

}

void cnn_visualization::write_scores(std::vector<float> scoreList, bool isRec)
{
    std::string filename;
    std::string molString;
    if(isRec)
    {
        filename = outRec;
        molString = recPDB;
    }
    else
    {
        filename = outLig;
        molString = ligPDB;
    }

    std::ofstream outFile; outFile.open(filename);

    outFile << "CNN MODEL: " << cnnopts.cnn_model << '\n';
    outFile << "CNN WEIGHTS: " << cnnopts.cnn_weights << '\n';

    std::stringstream molStream(molString);
    std::string line;
    std::string indexString;
    int index;
    std::stringstream scoreStream;
    std::string scoreString;
    while(std::getline(molStream, line))
    {
        if ((line.find("ATOM") < std::string::npos) || 
            (line.find("HETATM") < std::string::npos))
        {
            scoreStream.str(""); //clear stream for next score
            indexString = line.substr(6,5);
            index = std::stoi(indexString);

            if ((scoreList[index] > 0.001) || (scoreList[index] < -0.001)) //ignore very small scores
            {
                scoreStream << std::fixed << std::setprecision(5) << scoreList[index];
                outFile << line.substr(0,61);
                scoreString = scoreStream.str();
                scoreString.resize(5);
                outFile.width(5);
                outFile.fill('.');
                outFile << std::right << scoreString;
                outFile << line.substr(66) << '\n';

            }
            else
            {
                outFile << line << '\n';
            }
        }
        else
        {
            outFile << line << '\n';
        }

    }
}


bool cnn_visualization::check_in_range(std::set<int> atomList)
{
    float x = cenCoords[0];
    float y = cenCoords[1];
    float z = cenCoords[2];

    float allowedDist = visopts.box_size / 2;
    int numAtoms = recMol.NumAtoms();

    OBAtom* atom;
    for( auto i = atomList.begin(); i != atomList.end(); ++i)
    {
        if (*i >= numAtoms)
        {
            return false;
        }

        atom = recMol.GetAtom(*i);
        if(atom->GetX() < x + allowedDist)
            if (atom->GetY() < y + allowedDist)
                if (atom->GetZ() < z + allowedDist)
                    if (atom->GetX() > x - allowedDist)
                        if (atom->GetY() > y - allowedDist)
                            if (atom->GetZ() > z - allowedDist)
                                return true;
    }

        return false;
}

void cnn_visualization::ligCenter()
{
    vector3 cen = hLigMol.Center(0);
    cenCoords[0] = cen.GetX();
    cenCoords[1] = cen.GetY();
    cenCoords[2] = cen.GetZ();
}

std::vector<float> cnn_visualization::transform(std::vector<float> inList)
{
    std::vector<float> outList (inList.size());
    float tempVal;
    for (int i = 0; i < inList.size(); ++i)
    {
        tempVal = inList[i];
        if(tempVal < 0)
        {
            tempVal = 0 - std::sqrt(std::abs(tempVal));
        }
        else
        {
            tempVal = std::sqrt(tempVal);
        }

        tempVal = tempVal * 100;
        std::cout << inList[i] << " : " << tempVal << '\n';

        outList[i] = tempVal;
    }

    for (int i = 0; i < outList.size(); ++i)
    {
        std::cout << outList[i] << '\n';
    }

    return outList;
}
void cnn_visualization::remove_residues()
{
    std::vector<float> scoreDict(hRecMol.NumAtoms() + 1, 0.00);
    std::string lastRes = "";
    std::string currRes;
    std::vector<bool> atomList (hRecMol.NumAtoms() + 1, false);
    std::set<std::string> resList;

    std::string molString = hRec;
    std::stringstream molStream(molString);
    std::string line;
    while(std::getline(molStream, line))
    {
        if((line.find("ATOM") < std::string::npos) ||
           (line.find("HETATM") < std::string::npos))
        {
            currRes = line.substr(23,4);
            if(line.substr(23,4) != lastRes)
            {
                resList.insert(currRes);
                lastRes = currRes;
            }
        }
    }

    for( auto i = resList.begin(); i != resList.end(); ++i)
    {
        molStream.clear();
        molStream.str(molString);
        while(std::getline(molStream, line))
        {
            if((line.find("ATOM") < std::string::npos) ||
               (line.find("HETATM") < std::string::npos))
            {
                if(line.substr(23,4) == *i)
                {
                    std::string indexString = line.substr(6,5);
                    int index = std::stoi(indexString);
                    atomList[index] = true;
                }
            }
        }
        float scoreVal = remove_and_score(atomList, true);

        
        
        for ( auto f : atomList)
        {
            if(f)
            {
            scoreDict[f] = scoreVal;
            }
        }
        

        atomList.clear();
        atomList = std::vector<bool>(hRecMol.NumAtoms() + 1, false);

        
    }

    write_scores(scoreDict, true);
}

void cnn_visualization::remove_each_atom()
{
    std::vector<float> scoreDict(hLigMol.NumAtoms());
    std::stringstream ss (hLig);
    std::string line;

    std::string indexString;
    int index;
    std::vector<bool> removeList (hLigMol.NumAtoms() + 1);
    float scoreVal;

    while(std::getline(ss, line))
    {
      
        if ((line.find("ATOM") < std::string::npos) ||
            (line.find("HETATM") < std::string::npos))
        {
            indexString = line.substr(6, 5);
            index = std::stoi(indexString);
            if (hLigMol.GetAtom(index)->GetAtomicNum() != 1)
            {
                removeList[index] = true;

                scoreVal = remove_and_score(removeList, false);
                removeList[index] = false;

                scoreDict[index] = scoreVal;
            }
        }

    }

    write_scores(scoreDict, false);


}



