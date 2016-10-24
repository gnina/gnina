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

ColoredMol::ColoredMol (const vis_options &visopts, const cnn_options &cnnopts, FlexInfo &finfo, tee &log, const vec &center )
    {
        ligName = visopts.ligName;
        recName = visopts.recName;
        cnnmodel = cnnopts.cnn_model;
        weights = cnnopts.cnn_weights;
        size = visopts.size;
        outRec = visopts.outRec;
        outLig = visopts.outLig;
        frags_only = visopts.frags_only;
        atoms_only = visopts.atoms_only;
        verbose = visopts.verbose;
        this->cnnopts = cnnopts;
        std::cout << cnnopts.cnn_model;
        this->finfo = &finfo;
        this->log = &log;
        this->center = &center;
    }

void ColoredMol::color()
{
    OBConversion conv;

    std::string ext = boost::filesystem::extension(ligName);
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
        std::cout << "File extension not supported: " << ligName << '\n';
        std::cout << "Please use .pdb or .pdbqt for ligand\n";
        exit(0);
    }

    conv.ReadFile(&ligMol, ligName);

    ext = boost::filesystem::extension(recName);
    if(ext.compare(".pdb") == 0)
    {
      conv.SetInFormat("PDB");
    }
    else
    {
        std::cout << "File extension not supported: " << recName << '\n';
        std::cout << "Please use .pdb for receptor\n";
        exit(0);
    }

    conv.ReadFile(&recMol, recName);

    addHydrogens();
    ligCenter();

    //removeResidues();
    removeEachAtom();
}

void ColoredMol::print()
{
    std::cout << "ligName: " << ligName << '\n';
    std::cout << "recName: " << recName << '\n';
    std::cout << "cnnmodel: " << cnnmodel << '\n';
    std::cout << "size: " << size << '\n';
    std::cout << "outRec: " << outRec << '\n';
    std::cout << "outLig: " << outLig << '\n';
    std::cout << "frags_only: " << frags_only << '\n';
    std::cout << "atoms_only: " << atoms_only << '\n';
    std::cout << "verbose: " << verbose << '\n';
}

float ColoredMol::removeAndScore(std::vector<bool> removeList, bool isRec)
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

    if(!(isRec)) //if ligand
    {
        OBAtom* atom;
        for(int i = 0;i < removeList.size(); ++i)
        {
            if (removeList[i]) //index is in removeList
            {
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
        //make set for inRange test
        std::set<int> removeSet;
        for (int i = 0; i < removeList.size(); ++i)
        {
            if (removeList[i])
            {
                removeSet.insert(i);
            }
        }

       if (!(inRange(removeSet)))
        {
            return 0.00;
        }
    }

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
}
void ColoredMol::addHydrogens()
{
    recMol.AddHydrogens();
    ligMol.AddHydrogens();

    OBConversion conv;
    conv.SetInFormat("PDB");

    conv.SetOutFormat("PDB");
    recPDB = conv.WriteString(&recMol); //store rec pdb to preserve residue info

    
    conv.SetOutFormat("PDBQT"); //use pdbqt to make passing to parse_pdbqt easier
    conv.AddOption("r",OBConversion::OUTOPTIONS);
    conv.AddOption("c",OBConversion::OUTOPTIONS);
    hLig = conv.WriteString(&ligMol);
    hRec = conv.WriteString(&recMol);

    conv.SetInFormat("PDBQT");
    conv.ReadString(&hRecMol, hRec);
    conv.ReadString(&hLigMol, hLig);
}

float ColoredMol::score(const std::string &molString, bool isRec)
{
  if( !isRec )
  {
        std::stringstream ligStream(molString);
       // std::cout << molString << '\n';
        std::stringstream recStream(hRec);

        model m = parse_receptor_pdbqt("", recStream);

        //MolGetter mols;
        //mols.create_init_model(recName, "", *finfo, *log);

        CNNScorer cnn_scorer(cnnopts, *center, m);
        //mols.setInputFile(ligName);

        m.append(parse_ligand_stream_pdbqt(molString, ligStream));
        //bool worked = mols.readMoleculeIntoModel(*lig);
        //std::cout << "Worked: " << worked << '\n';

        float theScore = cnn_scorer.score(m);
        std::cout << "SCORE: " << theScore << '\n';

        return theScore;
        //return 1;
  }
  else
  {
    std::cout << "receptor not implemented, exiting\n";
    exit(0);
  }

}

void ColoredMol::writeScores(std::vector<float> scoreList, bool isRec)
{
    std::string filename;
    std::string molString;
    if(isRec)
    {
        filename = outRec;
        molString = hRec;
    }
    else
    {
        filename = outLig;
        molString = hLig;
    }

    std::ofstream outFile; outFile.open(filename);

    outFile << "CNN MODEL: " << cnnmodel << '\n';
    outFile << "CNN WEIGHTS: " << weights << '\n';

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


bool ColoredMol::inRange(std::set<int> atomList)
{
    float x = cenCoords[0];
    float y = cenCoords[1];
    float z = cenCoords[2];

    float allowedDist = size / 2;
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

void ColoredMol::ligCenter()
{
    vector3 cen = hLigMol.Center(0);
    cenCoords[0] = cen.GetX();
    cenCoords[1] = cen.GetY();
    cenCoords[2] = cen.GetZ();
}

std::vector<float> ColoredMol::transform(std::vector<float> inList)
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
void ColoredMol::removeResidues()
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
        float scoreVal = removeAndScore(atomList, true);

        
        
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

    writeScores(scoreDict, true);
}

void ColoredMol::removeEachAtom()
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

                scoreVal = removeAndScore(removeList, false);
                removeList[index] = false;

                scoreDict[index] = scoreVal;
            }
        }

    }


}



