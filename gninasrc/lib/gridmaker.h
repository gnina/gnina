/*
 * gridmaker.h
 *
 * Fills in grids with atom type information.  Implemented within the header
 * to make it easier to include as a dependency (cuda code separate though).
 *  Created on: May 31, 2016
 *      Author: dkoes
 */

#ifndef _GRIDMAKER_H_
#define _GRIDMAKER_H_

#include <vector>
#include <cmath>
#include <cuda.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <vector_types.h>
#include <boost/array.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/algorithm/string.hpp>

#ifndef VINA_ATOM_CONSTANTS_H
#include "gninasrc/lib/atom_constants.h"
#endif
#include "gpu_math.h"

using namespace std;

class GridMaker {
  boost::array< pair<float, float>, 3> dims;
  float3 center;
  float radiusmultiple;
  float resolution;
  float dimension;
  unsigned dim; //number of points on each side (cube)
	bool binary;

public:
	typedef boost::math::quaternion<float> quaternion;


	GridMaker(float res=0, float d=0, float rm = 1.5, bool b = false):
		radiusmultiple(rm), resolution(res), dimension(d), binary(b)
	{
	  initialize(res,d,rm,b);
	}

	virtual ~GridMaker() {}

	void initialize(float res, float d, float rm = 1.5, bool b = false) {
		resolution = res;
		dimension = d;
		radiusmultiple = rm;
		binary = b;
	  dim = ::round(dimension/resolution)+1; //number of grid points on a size
	  center.x = center.y = center.z = 0;
	}
	//mus set center before gridding
	void setCenter(double x, double y, double z)
	{
    center.x = x;
    center.y = y;
    center.z = z;
    float half = dimension/2.0;
    dims[0].first = x - half;
    dims[0].second = x + half;
    dims[1].first = y - half;
    dims[1].second = y + half;
    dims[2].first = z - half;
    dims[2].second = z + half;
	}

	template<typename Grid>
	void zeroGridsCPU(vector<Grid>& grids)
	{
		for (unsigned i = 0, n = grids.size(); i < n; i++)
		{
			std::fill(grids[i].data(), grids[i].data() + grids[i].num_elements(), 0.0);
		}
	}
  
  template<typename Grids>
  void zeroGridsCPU(Grids& grids)
  {
      std::fill(grids.data(), grids.data() + grids.num_elements(), 0.0);
  }
  
	pair<unsigned, unsigned> getrange(const pair<float, float>& d, double c, double r)
	{
	  pair<unsigned, unsigned> ret(0, 0);
	  double low = c - r - d.first;
	  if (low > 0)
	  {
	    ret.first = floor(low / resolution);
	  }

	  double high = c + r - d.first;
	  if (high > 0) //otherwise zero
	  {
	    ret.second = std::min(dim, (unsigned) ceil(high / resolution));
	  }
	  return ret;
	}

	//return the occupancy for atom a at point x,y,z
	float calcPoint(const float3& coords, double ar, float x, float y, float z)
	{
		float dx = x-coords.x;
		float dy = y-coords.y;
		float dz = z-coords.z;

	  float rsq = dx*dx+dy*dy+dz*dz;
	  if (binary)
	  {
	    //is point within radius?
	    if (rsq < ar * ar)
	      return 1.0;
	    else
	      return 0.0;
	  }
	  else
	  {
	    //for non binary we want a gaussian were 2 std occurs at the radius
	    //after which which switch to a quadratic
	    //the quadratic is to fit to have both the same value and first order
	    //derivative at the cross over point and a value and derivative of zero
	    //at 1.5*radius
	    double dist = sqrt(rsq);
	    if (dist >= ar * radiusmultiple)
	    {
	      return 0.0;
	    }
	    else if (dist <= ar)
	    {
	      //return gaussian
	      float h = 0.5 * ar;
	      float ex = -dist * dist / (2 * h * h);
	      return exp(ex);
	    }
	    else //return quadratic
	    {
	      float h = 0.5 * ar;
	      float eval = 1.0 / (M_E * M_E); //e^(-2)
	      float q = dist * dist * eval / (h * h) - 6.0 * eval * dist / h
	          + 9.0 * eval;
	      return q;
	    }
	  }
      //this was here, not sure why, giving "statement unreachable" compiler
      //warning and you know how much I love placating the compiler
	  //return 0.0;
	}


	//Grids is either a vector of multi_arrays or a multi_array
	//set the relevant grid points for passed info
	template<typename Grids>
	void setAtomsCPU(const vector<float4>& ainfo, const vector<short>& gridindex,  const quaternion& Q, Grids& grids)
	{
		zeroGridsCPU(grids);
		for (unsigned i = 0, n = ainfo.size(); i < n; i++)
		{
			int pos = gridindex[i];
			if (pos >= 0)
				setAtomCPU(ainfo[i], pos, Q, grids);
		}
	}

	//set the relevant grid points for provided atom
	template<typename Grids>
	void setAtomCPU(float4 ainfo, int whichgrid, const quaternion& Q, Grids& grids)
	{
		float radius = ainfo.w;
	  float r = radius * radiusmultiple;
	  float3 coords;
	  if (Q.real() != 0)
	  { //apply rotation
	    quaternion p(0, ainfo.x-center.x, ainfo.y-center.y, ainfo.z-center.z);
	    p = Q * p * (conj(Q) / norm(Q));

	    coords.x = p.R_component_2() + center.x;
	    coords.y = p.R_component_3() + center.y;
	    coords.z = p.R_component_4() + center.z;
	  }
	  else
	  {
	    coords.x = ainfo.x;
	    coords.y = ainfo.y;
	    coords.z = ainfo.z;
	  }

	  vector<pair<unsigned, unsigned> > ranges(3);
	  ranges[0] = getrange(dims[0], coords.x, r);
	  ranges[1] = getrange(dims[1], coords.y, r);
	  ranges[2] = getrange(dims[2], coords.z, r);

	  //for every grid point possibly overlapped by this atom
	  for (unsigned i = ranges[0].first, iend = ranges[0].second; i < iend; i++)
	  {
	    for (unsigned j = ranges[1].first, jend = ranges[1].second; j < jend;
	        j++)
	    {
	      for (unsigned k = ranges[2].first, kend = ranges[2].second;
	          k < kend; k++)
	      {
	        float x = dims[0].first + i * resolution;
	        float y = dims[1].first + j * resolution;
	        float z = dims[2].first + k * resolution;
	        float val = calcPoint(coords, radius, x, y, z);

	        if (binary)
	        {
	          if (val != 0)
	            grids[whichgrid][i][j][k] = 1.0; //don't add, just 1 or 0
	        }
	        else
	          grids[whichgrid][i][j][k] += val;

	      }
	    }
	  }
	}


	//GPU accelerated version, defined in cu file
	//pointers must point to GPU memory
	template<typename Dtype>
	void setAtomsGPU(unsigned natoms, float4 *coords, short *gridindex, quaternion Q, unsigned ngrids, Dtype *grids);


  void zeroAtomGradientsCPU(vector<float3>& agrad)
  {
    for (unsigned i = 0, n = agrad.size(); i < n; ++i)
    { 
      agrad[i].x = 0.0;
      agrad[i].y = 0.0;
      agrad[i].z = 0.0;
    }
  }

  //accumulate gradient from grid point x,y,z for provided atom
  void accumulateAtomGradient(const float3& coords, double ar, float x, float y, float z,
                              float gridval, float3& agrad)
  {
    float dist_x = x-coords.x;
    float dist_y = y-coords.y;
    float dist_z = z-coords.z;
    float dist2 = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
    if (dist2 == 0.0 || gridval == 0.0) // gradient is zero
      return;

    //in order to compute the loss gradient wrt atom position,
    //  we sum gradient from each gridpoint overlapped by the atom
    //  we can compute this as d_gridpoint/d_atompos * d_loss/d_gridpoint

    //the first term is a derivative of the atom-gridding operation:
    //  how does the occupancy value at grid position x,y,z change
    //  with respect to the coordinates of this atom?
    //the second term is gridval, computed by caffe backward pass to input layer

    //atom-gridding is circularly symmetric, so we compute the derivative
    //  relative to the distance from atom to gridpoint: d_gridpoint/d_atomdist
    //  then to convert to x,y,z, multiply by: d_atomdist/d_atomx (same for y,z)

    double dist = sqrt(dist2);
    float agrad_dist = 0.0;
    if (dist >= ar * radiusmultiple)
    { 
      return;
    }
    else if (dist <= ar) //gaussian derivative
    { 
      float h = 0.5 * ar;
      float ex = -dist2 / (2 * h * h);
      float coef = -dist / (h * h);
      agrad_dist = coef * exp(ex);
    }
    else //quadratic derivative
    { 
      float h = 0.5 * ar;
      float inv_e2 = 1.0 / (M_E * M_E); //e^(-2)
      agrad_dist = 2.0 * dist * inv_e2 / (h * h) - 6.0 * inv_e2 / h;
    }
    // d_loss/d_atomx = d_atomdist/d_atomx * d_gridpoint/d_atomdist * d_loss/d_gridpoint
    // sum across all gridpoints
    agrad.x += (-dist_x / dist) * agrad_dist * gridval;
    agrad.y += (-dist_y / dist) * agrad_dist * gridval;
    agrad.z += (-dist_z / dist) * agrad_dist * gridval;
  }

  //get the atom position gradient from relevant grid points for provided atom
  template<typename Grids>
  void setAtomGradientCPU(const float4& ainfo, int whichgrid, const quaternion& Q, const Grids& grids,
                          float3& agrad)
  { 
    float3 coords;
    if (Q.real() != 0) //apply rotation
    {
      quaternion p(0, ainfo.x-center.x, ainfo.y-center.y, ainfo.z-center.z);
      p = Q * p * (conj(Q) / norm(Q));

      coords.x = p.R_component_2() + center.x;
      coords.y = p.R_component_3() + center.y;
      coords.z = p.R_component_4() + center.z;
    }
    else
    {
      coords.x = ainfo.x;
      coords.y = ainfo.y;
      coords.z = ainfo.z;
    }

    //get grid index ranges that could possibly be overlapped by atom
    float radius = ainfo.w;
    float r = radius * radiusmultiple;
    vector<pair<unsigned, unsigned> > ranges(3);
    ranges[0] = getrange(dims[0], coords.x, r);
    ranges[1] = getrange(dims[1], coords.y, r);
    ranges[2] = getrange(dims[2], coords.z, r);

    //for every grid point possibly overlapped by this atom
    for (unsigned i = ranges[0].first, iend = ranges[0].second; i < iend; ++i)
    {
      for (unsigned j = ranges[1].first, jend = ranges[1].second; j < jend; ++j)
      {
        for (unsigned k = ranges[2].first, kend = ranges[2].second; k < kend; ++k)
        {
          //convert grid point coordinates to angstroms
          float x = dims[0].first + i * resolution;
          float y = dims[1].first + j * resolution;
          float z = dims[2].first + k * resolution;
          accumulateAtomGradient(coords, radius, x, y, z, grids[whichgrid][i][j][k], agrad);
        }
      }
    }
  }

  //backpropagate the gradient from atom grid to atom x,y,z positions
  template<typename Grids>
  void setAtomGradientsCPU(const vector<float4>& ainfo, const vector<short>& gridindex, const quaternion& Q,
                           const Grids& grids, vector<float3>& agrad)
  { 
    zeroAtomGradientsCPU(agrad);
    for (unsigned i = 0, n = ainfo.size(); i < n; ++i)
    {
      int whichgrid = gridindex[i]; // this is which atom-type channel of the grid to look at
      if (whichgrid >= 0)
        setAtomGradientCPU(ainfo[i], whichgrid, Q, grids, agrad[i]);
    }
  }


	static unsigned createDefaultMap(const char *names[], vector<int>& map)
	{
		map.assign(smina_atom_type::NumTypes, -1);
		const char **nameptr = names;
		unsigned cnt = 0;
		while (*nameptr != NULL)
		{
			string name(*nameptr);
			//note that if we every start using merged atom types by default
			//this code will have to be updated
			smt t = string_to_smina_type(name);
			if (t < smina_atom_type::NumTypes) //valid
			{
				map[t] = cnt;
				cnt++;
			}
			else //should never happen
			{
				cerr << "Invalid atom type " << name << "\n";
				exit(-1);
			}
			nameptr++;
		}
		return cnt;
	}

	//initialize default receptor/ligand maps
	//these were determined by an analysis of type frequencies
	static unsigned createDefaultRecMap(vector<int>& map)
	{
		const char *names[] =
		{ "AliphaticCarbonXSHydrophobe",
				"AliphaticCarbonXSNonHydrophobe",
				"AromaticCarbonXSHydrophobe",
				"AromaticCarbonXSNonHydrophobe",
				"Calcium",
				"Iron",
				"Magnesium",
				"Nitrogen",
				"NitrogenXSAcceptor",
				"NitrogenXSDonor",
				"NitrogenXSDonorAcceptor",
				"OxygenXSAcceptor",
				"OxygenXSDonorAcceptor",
				"Phosphorus",
				"Sulfur",
				"Zinc", NULL };

		return createDefaultMap(names, map);
	}

	static unsigned createDefaultLigMap(vector<int>& map)
	{
		const char *names[] =
		{ "AliphaticCarbonXSHydrophobe",
				"AliphaticCarbonXSNonHydrophobe",
				"AromaticCarbonXSHydrophobe",
				"AromaticCarbonXSNonHydrophobe",
				"Bromine",
				"Chlorine",
				"Fluorine",
				"Nitrogen",
				"NitrogenXSAcceptor",
				"NitrogenXSDonor",
				"NitrogenXSDonorAcceptor",
				"Oxygen",
				"OxygenXSAcceptor",
				"OxygenXSDonorAcceptor",
				"Phosphorus",
				"Sulfur",
				"SulfurAcceptor",
				"Iodine",
				NULL };
		return createDefaultMap(names, map);
	}

	//create a mapping from atom type ids to a unique id given a file specifying
	//what types we care about (anything missing is ignored); if multiple types are
	//on the same line, they are merged, if the file isn't specified, use default mapping
	//return total number of types
	//map is indexed by smina_atom_type, maps to -1 if type should be ignored
	static unsigned createAtomTypeMap(const string& fname, vector<int>& map)
	{
		using namespace std;
		using namespace boost::algorithm;
	  map.assign(smina_atom_type::NumTypes, -1);

	  if(fname.size() == 0) {
	  	std::cerr <<  "Map file not specified\n";
	  	exit(-1);
	  }

		unsigned cnt = 0;
		ifstream in(fname.c_str());

		if(!in) {
			std::cerr << "Could not open " << fname << "\n";
			exit(-1);
		}

		string line;
		while (getline(in, line))
		{
			vector<string> types;
			split(types, line, is_any_of("\t \n"));
			for (unsigned i = 0, n = types.size(); i < n; i++)
			{
				const string& name = types[i];
				smt t = string_to_smina_type(name);
				if (t < smina_atom_type::NumTypes) //valid
				{
					map[t] = cnt;
				}
				else if (name.size() > 0) //this ignores consecutive delimiters
				{
					cerr << "Invalid atom type " << name << "\n";
					exit(-1);
				}
			}
			if (types.size() > 0)
				cnt++;
		}
		return cnt;
	}

};

#endif /* _GRIDMAKER_H_ */
