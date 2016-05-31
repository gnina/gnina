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
	typedef boost::math::quaternion<double> quaternion;


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
	  dim = round(dimension/resolution)+1; //number of grid points on a size
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

	template<typename Grids>
	void zeroGridsCPU(Grids& grids)
	{
		for (unsigned i = 0, n = grids.size(); i < n; i++)
		{
			std::fill(grids[i].data(), grids[i].data() + grids[i].num_elements(), 0.0);
		}
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
	  return 0.0;
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


	//GPU accelerated version
	void setAtomsGPU(unsigned natoms, float4 *coords, short *gridindex, unsigned ngrids, float *grids);
};

#endif /* _GRIDMAKER_H_ */
