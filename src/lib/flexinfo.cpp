#include "flexinfo.h"
#include <boost/lexical_cast.hpp>
#include <openbabel/mol.h>
#include <boost/unordered_map.hpp>

FlexInfo::FlexInfo(const std::string& flexres, double flexdist, const std::string& ligand, tee& l): flex_dist(flexdist), log(l)
{
	using namespace OpenBabel;
	using namespace std;
	//first extract comma separated list
	if(flexres.size() > 0)
	{
		vector<string> tokens;
		boost::split(tokens, flexres, boost::is_any_of(","));

		vector<string> chres;
		for(unsigned i = 0, n = tokens.size(); i < n; i++)
		{
			//each token may be either chain:resid or just resid
			string tok = tokens[i];
			boost::split(chres, tok, boost::is_any_of(":"));
			char chain = 0;
			int resid = 0;
			if(chres.size() == 2)
			{
				if(chres[0].size() != 1)
				log << "WARNING: chain specification not single character " << chres[0] << "\n";
				chain = chres[0][0]; //if empty will be null which is what is desired
				resid = boost::lexical_cast<int>(chres[1]);
			}
			else if(chres.size() == 1)
			{
				resid = boost::lexical_cast<int>(chres[0]);
			}
			else
			{
				log << "WARNING: ignoring invalid chain:resid specifier " << tok << "\n";
				continue;
			}

			residues.insert(pair<char,int>(chain,resid));
		}
	}

	if(ligand.size() > 0 && flex_dist > 0)
	{
		//next read ligand for distance checking
		obmol_opener opener;
		OBConversion conv;
		opener.openForInput(conv, ligand);
		conv.Read(&distligand);//first ligand only
	}

}

void FlexInfo::extractFlex(OpenBabel::OBMol& receptor, OpenBabel::OBMol& rigid,
		std::string& flexpdbqt)
{
	using namespace OpenBabel;
	rigid = receptor;

	flexpdbqt.clear();

	//identify residues close to distligand here
	Box b;
	b.add_ligand_box(distligand);
	b.expand(flex_dist);
	double flsq = flex_dist * flex_dist;

	FOR_ATOMS_OF_MOL(a, rigid)
	{
		if(a->GetAtomicNum() == 1)
			continue; //heavy atoms only
		vector3 v = a->GetVector();
		if (b.ptIn(v.x(), v.y(), v.z()))
		{
			//in box, see if any atoms are close enough
			FOR_ATOMS_OF_MOL(b, distligand)
			{
				vector3 bv = b->GetVector();
				if (v.distSq(bv) < flsq)
				{
					//process residue
					OBResidue *residue = a->GetResidue();
					if (residue)
					{
						char ch = residue->GetChain();
						int resid = residue->GetNum();
						residues.insert(std::pair<char, int>(ch, resid));
					}
					break;
				}
			}
		}
	}

	//replace any empty chains with first chain in mol
	char defaultch = ' ';
	OBResidueIterator ritr;
	OBResidue *firstres = rigid.BeginResidue(ritr);
	if (firstres)
		defaultch = firstres->GetChain();

	std::vector<std::pair<char, int> > sortedres(residues.begin(),
			residues.end());
	for (unsigned i = 0, n = sortedres.size(); i < n; i++)
	{
		if (sortedres[i].first == 0)
			sortedres[i].first = defaultch;
	}

	sort(sortedres.begin(), sortedres.end());

	if (sortedres.size() > 0)
	{
		log << "Flexible residues:";
		for (unsigned i = 0, n = sortedres.size(); i < n; i++)
		{
			log << " " << sortedres[i].first << ":" << sortedres[i].second;
		}
		log << "\n";
	}
	//reinsert residues now with defaul chain
	residues.clear();
	residues.insert(sortedres.begin(), sortedres.end());

	OBConversion conv;
	conv.SetOutFormat("PDBQT");
	conv.AddOption("s", OBConversion::OUTOPTIONS); //flexible residue
	//identify atoms that have to be in flexible component
	//this is the side chain and CA, but _not_ the C and N
	for(OBResidueIterator ritr = rigid.BeginResidues(), rend = rigid.EndResidues(); ritr != rend; ++ritr)
	{
		OBResidue *r = *ritr;
		char ch = r->GetChain();
		int resid = r->GetNum();
		std::pair<char,int> chres(ch,resid);
		if(residues.count(chres))
		{
			//create a separate molecule for each flexible residue
			OBMol flex;
			std::vector<OBAtom*> flexatoms; //rigid atom ptrs that should be flexible
			boost::unordered_map<OBAtom*, int> flexmap; //map rigid atom ptrs to atom indices in flex

			//make absolutely sure that CA is the first atom
			for(OBAtomIterator aitr = r->BeginAtoms(), aend = r->EndAtoms(); aitr != aend; ++aitr)
			{
				OBAtom *a = *aitr;
				std::string aid = r->GetAtomID(a);
				boost::trim(aid);
				if(aid == "CA")
				{
					flexatoms.push_back(a);
					flex.AddAtom(*a);
					flexmap[a] = flex.NumAtoms(); //after addatom since indexed by
				}
			}

			for(OBAtomIterator aitr = r->BeginAtoms(), aend = r->EndAtoms(); aitr != aend; ++aitr)
			{
				OBAtom *a = *aitr;
				std::string aid = r->GetAtomID(a);
				boost::trim(aid);
				if(aid != "CA" && aid != "N" && aid != "C" &&
						aid != "O" && aid != "H" && //leave backbone alone other than CA which is handled above
						 !a->IsNonPolarHydrogen())
				{
					flexatoms.push_back(a);
					flex.AddAtom(*a);
					flexmap[a] = flex.NumAtoms(); //after addatom since indexed by
				}
			}

			//now add bonds - at some point I should add functions to openbabel to do this..
			for(unsigned i = 0, n = flexatoms.size(); i < n; i++)
			{
				OBAtom *a = flexatoms[i];
				FOR_BONDS_OF_ATOM(b, a)
				{
					OBBond& bond = *b;
					//just do one direction, if first atom is a
					//and second atom is a flexatom need to add bond
					if(a == bond.GetBeginAtom() && flexmap.count(bond.GetEndAtom()))
					{
						flex.AddBond(flexmap[a],flexmap[bond.GetEndAtom()],bond.GetBO(), bond.GetFlags());
					}
				}
			}

			flex.AddResidue(*r);
			OBResidue *newres = flex.GetResidue(0);
			if(newres)
			{
				//add all atoms with proper atom ids
				for(unsigned i = 0, n = flexatoms.size(); i < n; i++)
				{
					OBAtom *origa = flexatoms[i];
					OBAtom *newa = flex.GetAtom(flexmap[origa]);
					newres->AddAtom(newa);
					newa->SetResidue(newres);
					std::string aid = r->GetAtomID(origa);
					newres->SetAtomID(newa, aid);
				}
			}
			flexpdbqt += conv.WriteString(&flex);

			//remove flexatoms from rigid
			for(unsigned i = 0, n = flexatoms.size(); i < n; i++)
			{
				OBAtom *a = flexatoms[i];
				rigid.DeleteAtom(a);
			}
		} //end if residue
	}

}
