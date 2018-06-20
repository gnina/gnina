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
#include <stdarg.h>
#include <device_types.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <vector_types.h>
#include <boost/array.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/algorithm/string.hpp>
#include "quaternion.h"
#include "gridoptions.h"
#include "caffe/proto/caffe.pb.h"

#ifndef VINA_ATOM_CONSTANTS_H
#include "gninasrc/lib/atom_constants.h"
#endif
#include "gpu_math.h"

using namespace std;

class GridMaker {
  protected:
    float2 dims[3];
    float3 center;
    float radiusmultiple;
    float resolution;
    float dimension;
    float rsq; //radius squared (dimension/2)^2
    unsigned dim; //number of points on each side (grid)
    bool binary;
    bool spherize; //mask out atoms not within sphere of center
  
  public:
    typedef boost::math::quaternion<float> quaternion;
  
  
    GridMaker(float res = 0, float d = 0, float rm = 1.5, bool b = false,
        bool s = false)
        : radiusmultiple(rm), resolution(res), dimension(d), binary(b),
            spherize(s) {
      initialize(res,d,rm,b,s);
    }
  
    virtual ~GridMaker() {}
  
    virtual void initialize(float res, float d, float rm = 1.5, bool b = false, bool s =
        false) {
      resolution = res;
      dimension = d;
      radiusmultiple = rm;
      binary = b;
      spherize = s;
      dim = ::round(dimension / resolution) + 1; //number of grid points on a side
      rsq = (dimension / 2) * (dimension / 2);
      center.x = center.y = center.z = 0;
    }
  
    virtual void initialize(const caffe::MolGridDataParameter& param) {
      initialize(param.resolution(), param.dimension(), param.radius_multiple(), 
          param.binary_occupancy(), param.spherical_mask());
    }
  
    virtual void initialize(const gridoptions& opt, float rm) {
      initialize(opt.res, opt.dim, rm, opt.binary, opt.spherize);
    }
  
    //mus set center before gridding
    virtual void setCenter(double x, double y, double z) {
      center.x = x;
      center.y = y;
      center.z = z;
      float half = dimension/2.0;
      dims[0].x = x - half;
      dims[0].y = x + half;
      dims[1].x = y - half;
      dims[1].y = y + half;
      dims[2].x = z - half;
      dims[2].y = z + half;
    }
  
    template<typename Grid>
    void zeroGridsCPU(vector<Grid>& grids) {
      for (unsigned i = 0, n = grids.size(); i < n; i++) {
        std::fill(grids[i].data(), grids[i].data() + grids[i].num_elements(), 0.0);
      }
    }
    
    template<typename Grids>
    void zeroGridsCPU(Grids& grids) {
      std::fill(grids.data(), grids.data() + grids.num_elements(), 0.0);
    }

    pair<unsigned, unsigned> getrange(const float2& d, double c, double r) {
      pair<unsigned, unsigned> ret(0, 0);
      double low = c - r - d.x;
      if (low > 0) {
        ret.first = floor(low / resolution);
      }

      double high = c + r - d.x;
      if (high > 0) //otherwise zero
          {
        ret.second = std::min(dim, (unsigned) ceil(high / resolution));
      }
      return ret;
    }

    __device__ uint2 getrange_gpu(const float2& d, double c, double r) {
      uint2 ret = make_uint2(0, 0);
      double low = c - r - d.x;

      if (low > 0) {
        ret.x = floor(low / resolution);
      }

      double high = c + r - d.x;
      if (high > 0) //otherwise zero
          {
        ret.y = min(dim, (unsigned) ceil(high / resolution));
      }
      return ret;
    }

    //return the occupancy for atom a at point x,y,z
    float calcPoint(const float3& coords, double ar, float x, float y,
        float z) {
      float dx = x - coords.x;
      float dy = y - coords.y;
      float dz = z - coords.z;

      float rsq = dx * dx + dy * dy + dz * dz;
      if (binary) {
        //is point within radius?
        if (rsq < ar * ar)
          return 1.0;
        else
          return 0.0;
      } else {
        //for non binary we want a gaussian were 2 std occurs at the radius
        //after which which switch to a quadratic
        //the quadratic is to fit to have both the same value and first order
        //derivative at the cross over point and a value and derivative of zero
        //at 1.5*radius
        double dist = sqrt(rsq);
        if (dist >= ar * radiusmultiple) {
          return 0.0;
        } else
          if (dist <= ar) {
            //return gaussian
            float h = 0.5 * ar;
            float ex = -dist * dist / (2 * h * h);
            return exp(ex);
          } else //return quadratic
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

    template <typename Dtype>
    void setAtomsCPU(const vector<float4>& ainfo, const vector<short>& gridindex,  
        const quaternion& Q, Dtype* data, unsigned ntypes) {
      boost::multi_array_ref<Dtype, 4> grids(data, boost::extents[ntypes][dim][dim][dim]);
      setAtomsCPU(ainfo, gridindex, Q, grids);
    }
  
    //Grids is either a vector of multi_arrays or a multi_array
    //set the relevant grid points for passed info
    template<typename Grids>
    void setAtomsCPU(const vector<float4>& ainfo, const vector<short>& gridindex,  
        const quaternion& Q, Grids& grids) {
      zeroGridsCPU(grids);
      for (unsigned i = 0, n = ainfo.size(); i < n; i++) {
        int pos = gridindex[i];
        if (pos >= 0)
          setAtomCPU(ainfo[i], pos, Q, grids);
      }
    }
  
    //set the relevant grid points for provided atom
    template<typename Grids>
    void setAtomCPU(float4 ainfo, int whichgrid, const quaternion& Q, 
        Grids& grids) {
      float radius = ainfo.w;
      float r = radius * radiusmultiple;
      float3 coords;
  
      if(spherize) {
        float xdiff = ainfo.x - center.x;
        float ydiff = ainfo.y - center.y;
        float zdiff = ainfo.z - center.z;
        float distsq = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
        if(distsq > rsq) {
          return; //ignore
        }
      }
  
      if (Q.real() != 0) { //apply rotation
        quaternion p(0, ainfo.x - center.x, ainfo.y - center.y, ainfo.z - center.z);
        p = Q * p * (conj(Q) / norm(Q));
        coords.x = p.R_component_2() + center.x;
        coords.y = p.R_component_3() + center.y;
        coords.z = p.R_component_4() + center.z;
      }
      else {
        coords.x = ainfo.x;
        coords.y = ainfo.y;
        coords.z = ainfo.z;
      }
  
      vector<pair<unsigned, unsigned> > ranges(3);
      ranges[0] = getrange(dims[0], coords.x, r);
      ranges[1] = getrange(dims[1], coords.y, r);
      ranges[2] = getrange(dims[2], coords.z, r);
  
      //for every grid point possibly overlapped by this atom
      for (unsigned i = ranges[0].first, iend = ranges[0].second; i < iend; i++) {
        for (unsigned j = ranges[1].first, jend = ranges[1].second; j < jend;
            j++) {
          for (unsigned k = ranges[2].first, kend = ranges[2].second;
              k < kend; k++) {
            float x = dims[0].x + i * resolution;
            float y = dims[1].x + j * resolution;
            float z = dims[2].x + k * resolution;
            float val = calcPoint(coords, radius, x, y, z);
  
            if (binary) {
              if (val != 0) 
                grids[whichgrid][i][j][k] = 1.0;
            }
            else {
              grids[whichgrid][i][j][k] += val;
            }
  
          }
        }
      }
    }
  
  
    //GPU accelerated version, defined in cu file
    //pointers must point to GPU memory
    template<typename Dtype>
    void setAtomsGPU(unsigned natoms, float4 *coords, short *gridindex, qt Q, 
        unsigned ngrids, Dtype *grids);
  
  
    void zeroAtomGradientsCPU(vector<float3>& agrad) {
      for (unsigned i = 0, n = agrad.size(); i < n; ++i) { 
        agrad[i].x = 0.0;
        agrad[i].y = 0.0;
        agrad[i].z = 0.0;
      }
    }
  
    __host__ __device__
    void accumulateAtomGradient(const float3& coords, double ar, float x,
              float y, float z, float gridval, float3& agrad, int whichgrid) {
      //sum gradient grid values overlapped by the atom times the
      //derivative of the atom density at each grid point
      float dist_x = x - coords.x;
      float dist_y = y - coords.y;
      float dist_z = z - coords.z;
      float dist2 = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
      double dist = sqrt(dist2);
      float agrad_dist = 0.0;
      if (dist >= ar * radiusmultiple) {//no overlap
        return;
      }
      else if (dist <= ar) {//gaussian derivative
        float h = 0.5 * ar;
        float ex = -dist2 / (2 * h * h);
        float coef = -dist / (h * h);
        agrad_dist = coef * exp(ex);
      }
      else {//quadratic derivative
        float h = 0.5 * ar;
        float inv_e2 = 1.0 / (M_E * M_E); //e^(-2)
        agrad_dist = 2.0 * dist * inv_e2 / (h * h) - 6.0 * inv_e2 / h;
      }
      // d_loss/d_atomx = d_atomdist/d_atomx * d_gridpoint/d_atomdist * d_loss/d_gridpoint
      // sum across all gridpoints
      //dkoes - the negative sign is because we are considering the derivative of the center vs grid
      float gx = -(dist_x / dist) * agrad_dist * gridval;
      float gy = -(dist_y / dist) * agrad_dist * gridval;
      float gz = -(dist_z / dist) * agrad_dist * gridval;
      agrad.x += gx;
      agrad.y += gy;
      agrad.z += gz;
    }
  
    //get the atom position gradient from relevant grid points for provided atom
    //if isrelevance is true, simply sum overlapping values
    template<typename Grids>
    void setAtomGradientCPU(const float4& ainfo, int whichgrid,
        const quaternion& Q, const Grids& grids,
        float3& agrad, bool isrelevance = false) {
      float3 coords;
      if (Q.real() != 0) {//apply rotation
        quaternion p(0, ainfo.x - center.x, ainfo.y - center.y,
            ainfo.z - center.z);
        p = Q * p * (conj(Q) / norm(Q));
  
        coords.x = p.R_component_2() + center.x;
        coords.y = p.R_component_3() + center.y;
        coords.z = p.R_component_4() + center.z;
      } 
      else {
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
      for (unsigned i = ranges[0].first, iend = ranges[0].second; i < iend;
          ++i) {
        for (unsigned j = ranges[1].first, jend = ranges[1].second;
            j < jend; ++j) {
          for (unsigned k = ranges[2].first, kend = ranges[2].second;
              k < kend; ++k) {
            //convert grid point coordinates to angstroms
            float x = dims[0].x + i * resolution;
            float y = dims[1].x + j * resolution;
            float z = dims[2].x + k * resolution;
  
            if (isrelevance) 
              accumulateAtomRelevance(coords, radius, x, y, z,
                  grids[whichgrid][i][j][k], agrad);
            else //true gradient, distance matters
              accumulateAtomGradient(coords, radius, x, y, z,
                  grids[whichgrid][i][j][k], agrad, whichgrid);
          }
        }
      }
    }
  
    template <typename Dtype>
    void setAtomGradientsCPU(const vector<float4>& ainfo, const vector<short>& gridindex, 
                             quaternion Q, Dtype* data, vector<float3>& agrad, 
                             unsigned offset, unsigned ntypes) {
      boost::multi_array_ref<Dtype, 4> grids(data+offset, boost::extents[ntypes][dim][dim][dim]);
      setAtomGradientsCPU(ainfo, gridindex, Q, grids, agrad);
    }
  
    //backpropagate the gradient from atom grid to atom x,y,z positions
    template<typename Grids>
    void setAtomGradientsCPU(const vector<float4>& ainfo, 
        const vector<short>& gridindex, const quaternion& Q, const Grids& grids, 
        vector<float3>& agrad) { 
      zeroAtomGradientsCPU(agrad);
      for (unsigned i = 0, n = ainfo.size(); i < n; ++i) {
        int whichgrid = gridindex[i]; // this is which atom-type channel of the grid to look at
        if (whichgrid >= 0) {
          setAtomGradientCPU(ainfo[i], whichgrid, Q, grids, agrad[i]);
        }
      }
    }
 
    __host__ __device__
    void accumulateAtomRelevance(const float3& coords, double ar, float x,
        float y, float z, float gridval, float3& agrad) {
      //simple sum of values that the atom overlaps
      float dist_x = x - coords.x;
      float dist_y = y - coords.y;
      float dist_z = z - coords.z;
      float dist2 = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
      double dist = sqrt(dist2);
      if (dist >= ar * radiusmultiple) {
        return;
      } else {
        agrad.x += gridval;
      }
    }

    template<typename Dtype>
    __device__
    void setAtomGradientsGPU(const float4* ainfo, short* gridindices, 
    float3* agrads, const qt Q, const Dtype* grids, unsigned remainder_offset, 
    bool isrelevance = false) {
#ifdef __CUDA_ARCH__
      //TODO: implement
      assert(isrelevance == false);

      int idx = blockDim.x * blockIdx.x + threadIdx.x + remainder_offset;
      int whichgrid = gridindices[idx];
      float4 atom = ainfo[idx];
      float3 coords; 

      if (Q.real() != 0) //apply rotation
          {
        float3 p = Q.rotate(atom.x - center.x, atom.y - center.y,
            atom.z - center.z);
        coords = p + center;
      } else {
        coords.x = atom.x;
        coords.y = atom.y;
        coords.z = atom.z;
      }

      //get grid index ranges that could possibly be overlapped by atom
      float radius = atom.w;
      float r = radius * radiusmultiple;
      uint2 ranges[3];
      ranges[0] = getrange_gpu(dims[0], coords.x, r);
      ranges[1] = getrange_gpu(dims[1], coords.y, r);
      ranges[2] = getrange_gpu(dims[2], coords.z, r);

      for (unsigned i = ranges[0].x, iend = ranges[0].y; i < iend; ++i) {
        for (unsigned j = ranges[1].x, jend = ranges[1].y; j < jend; ++j) {
          for (unsigned k = ranges[2].x, kend = ranges[2].y; k < kend; ++k) {
            //convert grid point coordinates to angstroms
            float x = dims[0].x + i * resolution;
            float y = dims[1].x + j * resolution;
            float z = dims[2].x + k * resolution;

            if (isrelevance) {
              accumulateAtomRelevance(coords, radius, x, y, z,
                  grids[(((whichgrid * dim) + i) * dim + j) * dim + k],
                  agrads[idx]);
            } else {
            accumulateAtomGradient(coords, radius, x, y, z, 
                      grids[(((whichgrid * dim) + i) * dim + j) * dim + k], 
                          agrads[idx], whichgrid);
            }
          }
        }
      }
#endif
    }

    //summ up gradient values overlapping atoms
    template<typename Grids>
    void setAtomRelevanceCPU(const vector<float4>& ainfo, 
        const vector<short>& gridindex, const quaternion& Q, const Grids& grids, 
        vector<float3>& agrad, unsigned batch_idx) {
      zeroAtomGradientsCPU(agrad);
      for (unsigned i = 0, n = ainfo.size(); i < n; ++i) {
        int whichgrid = gridindex[i]; // this is which atom-type channel of the grid to look at
        if (whichgrid >= 0) {
          setAtomGradientCPU(ainfo[i], whichgrid, Q, grids, agrad[i], true);
        }
      }
    }

    virtual unsigned createDefaultMap(const char *names[], vector<int>& map) {
      map.assign(smina_atom_type::NumTypes, -1);
      const char **nameptr = names;
      unsigned cnt = 0;
      while (*nameptr != NULL) {
        string line(*nameptr);
        vector<string> names;
        boost::algorithm::split(names, line, boost::is_space(), 
            boost::algorithm::token_compress_on);
        for(unsigned i = 0, n = names.size(); i < n; i++) {
          string name = names[i];
          smt t = string_to_smina_type(name);
          if(t < smina_atom_type::NumTypes) //valid
          {
            map[t] = cnt;
          }
          else //should never happen
          {
            cerr << "Invalid atom type " << name << "\n";
            exit(-1);
          }
        }

        if(names.size()) //skip empty lines
          cnt++;

        nameptr++;
      }
      return cnt;
    }

    //initialize default receptor/ligand maps
    //these were determined by an analysis of type frequencies
    virtual unsigned createDefaultRecMap(vector<int>& map) {
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

    virtual unsigned createDefaultLigMap(vector<int>& map) {
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
          "Boron",
          NULL };
      return createDefaultMap(names, map);
    }

    //create a mapping from atom type ids to a unique id given a file specifying
    //what types we care about (anything missing is ignored); if multiple types are
    //on the same line, they are merged, if the file isn't specified, use default mapping
    //return total number of types
    //map is indexed by smina_atom_type, maps to -1 if type should be ignored
    virtual unsigned createAtomTypeMap(const string& fname, vector<int>& map) {
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
      while (getline(in, line)) {
        vector<string> types;
        split(types, line, is_any_of("\t \n"));
        for (unsigned i = 0, n = types.size(); i < n; i++) {
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

class RNNGridMaker : public GridMaker {
  public:
    float subgrid_dim;
    unsigned batch_size;
    unsigned batch_idx;
    unsigned ntypes;
    unsigned nrec_types;
    unsigned nlig_types;
    RNNGridMaker(float res=0, float d=0, float rm=1.5, bool b=false, 
        bool s=false, float sd=3.0, unsigned bs=1, unsigned bi=0, 
        unsigned nt=0, unsigned nrt=0, unsigned nlt=0) : 
      GridMaker(res, d, rm, b, s), subgrid_dim(sd), batch_size(bs), 
      batch_idx(bi), ntypes(nt), nrec_types(nrt), nlig_types(nlt) {
      initialize(res, d, rm, b, s, sd, bs, bi, nt, nrt, nlt);
    }

    virtual ~RNNGridMaker() {}

    virtual void initialize(float res, float d, float rm=1.5, bool b = false, 
        bool s = false, float sd=3.0, unsigned bs=1, unsigned bi=0, 
        unsigned nt=0, unsigned nrt=0, unsigned nlt=0) {
      subgrid_dim = sd;
      batch_size = bs;
      batch_idx = bi;
      ntypes = nt;
      nrec_types = nrt;
      nlig_types = nlt;
      GridMaker::initialize(res, d, rm, b, s);
      if (subgrid_dim && fmod((dimension - subgrid_dim), subgrid_dim + 
            resolution) != 0)
        printf("Subgrid dimension must evenly divide total grid dimension");
    }

    virtual void initialize(const caffe::MolGridDataParameter& param) {
      initialize(param.resolution(), param.dimension(), param.radius_multiple(), 
          param.binary_occupancy(), param.spherical_mask(), param.subgrid_dim(), 
          param.batch_size());
    }

    virtual void initialize(const gridoptions& opt, float rm) {
      initialize(opt.res, opt.dim, rm, opt.binary, opt.spherize, opt.subgrid_dim);
    }

    virtual void setCenter(double x, double y, double z) {
      batch_idx = (batch_idx + 1) % batch_size;
      GridMaker::setCenter(x, y, z);
    }

    virtual unsigned createDefaultMap(const char *names[], vector<int>& map) {
      unsigned _ntypes = GridMaker::createDefaultMap(names, map);
      ntypes += _ntypes;
      return _ntypes;
    }

    virtual unsigned createDefaultRecMap(vector<int>& map) {
      nrec_types = GridMaker::createDefaultRecMap(map);
      return nrec_types;
    }

    virtual unsigned createDefaultLigMap(vector<int>& map) {
      nlig_types = GridMaker::createDefaultLigMap(map);
      return nlig_types;
    }

    virtual unsigned createAtomTypeMap(const string& fname, vector<int>& map) {
      unsigned _ntypes = GridMaker::createAtomTypeMap(fname, map);
      ntypes += _ntypes;
      return _ntypes;
    }

    template <typename Dtype>
    void setAtomsCPU(const vector<float4>& ainfo, const vector<short>& gridindex, 
        const quaternion& Q, Dtype* data, unsigned ntypes) {
      unsigned grids_per_dim = (this->dimension - subgrid_dim) / (subgrid_dim 
          + resolution) + 1;
      unsigned ngrids = grids_per_dim * grids_per_dim * grids_per_dim;
      unsigned subgrid_dim_in_points = dim / grids_per_dim;
      boost::multi_array_ref<Dtype, 6> grids(data, 
          boost::extents[ngrids][batch_size][ntypes][subgrid_dim_in_points]
          [subgrid_dim_in_points][subgrid_dim_in_points]);
      zeroGridsCPU(grids);
      for (unsigned i = 0, n = ainfo.size(); i < n; i++) {
        int pos = gridindex[i];
        if (pos >= 0)
          setAtomCPU(ainfo[i], pos, Q, grids);
      }
    }

    template<typename Grids>
    void setAtomsCPU(const vector<float4>& ainfo, 
        const vector<short>& gridindex,  const quaternion& Q, Grids& grids) {
      zeroGridsCPU(grids);
      for (unsigned i = 0, n = ainfo.size(); i < n; i++) {
        int pos = gridindex[i];
        if (pos >= 0)
          setAtomCPU(ainfo[i], pos, Q, grids);
      }
    }

    template<typename Grids>
    auto getGridElement(Grids& grids, unsigned grid_idx, unsigned whichgrid, 
        unsigned x, unsigned y, unsigned z) -> 
        decltype(&grids[0][0][0][0][0][0]){
      return &grids[grid_idx][batch_idx][whichgrid][x][y][z];
    }

    template<typename Allocator, typename Dtype>
    Dtype* getGridElement(std::vector<boost::multi_array<Dtype, 3, 
        Allocator>>& grids, unsigned grid_idx, unsigned whichgrid, unsigned x, 
        unsigned y, unsigned z) {
      return &grids[grid_idx * ntypes + whichgrid][x][y][z];
    }

    //TODO: possible to merge this with base version?
    template<typename Grids>
    void setAtomCPU(float4 ainfo, int whichgrid, const quaternion& Q, 
        Grids& grids) {
      float radius = ainfo.w;
      float r = radius * radiusmultiple;
      float3 coords;

      if(spherize) {
        float xdiff = ainfo.x - center.x;
        float ydiff = ainfo.y - center.y;
        float zdiff = ainfo.z - center.z;
        float distsq = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
        if(distsq > rsq) {
          return; //ignore
        }
      }

      if (Q.real() != 0) { //apply rotation
        quaternion p(0, ainfo.x - center.x, ainfo.y - center.y, ainfo.z - 
            center.z);
        p = Q * p * (conj(Q) / norm(Q));
        coords.x = p.R_component_2() + center.x;
        coords.y = p.R_component_3() + center.y;
        coords.z = p.R_component_4() + center.z;
      }
      else {
        coords.x = ainfo.x;
        coords.y = ainfo.y;
        coords.z = ainfo.z;
      }

      vector<pair<unsigned, unsigned> > ranges(3);
      ranges[0] = getrange(dims[0], coords.x, r);
      ranges[1] = getrange(dims[1], coords.y, r);
      ranges[2] = getrange(dims[2], coords.z, r);

      //for every grid point possibly overlapped by this atom
      for (unsigned i = ranges[0].first, iend = ranges[0].second; 
          i < iend; i++) {
        for (unsigned j = ranges[1].first, jend = ranges[1].second; 
            j < jend; j++) {
          for (unsigned k = ranges[2].first, kend = ranges[2].second;
              k < kend; k++) {
            float x = dims[0].x + i * resolution;
            float y = dims[1].x + j * resolution;
            float z = dims[2].x + k * resolution;
            float val = calcPoint(coords, radius, x, y, z);
            unsigned grids_per_dim = (dimension - subgrid_dim) / 
              (subgrid_dim + resolution) + 1;
            unsigned subgrid_dim_in_points = dim / grids_per_dim;
            unsigned subgrid_idx_x = i / subgrid_dim_in_points; 
            unsigned subgrid_idx_y = j / subgrid_dim_in_points; 
            unsigned subgrid_idx_z = k / subgrid_dim_in_points; 
            unsigned rel_x = i % subgrid_dim_in_points; 
            unsigned rel_y = j % subgrid_dim_in_points; 
            unsigned rel_z = k % subgrid_dim_in_points; 
            unsigned grid_idx = (((subgrid_idx_x * grids_per_dim) + 
                  subgrid_idx_y) * grids_per_dim + subgrid_idx_z);

            if (binary) {
              if (val != 0)
                *getGridElement(grids, grid_idx, whichgrid, rel_x, rel_y, 
                    rel_z) = 1.0;
            }
            else {
              *getGridElement(grids, grid_idx, whichgrid, rel_x, rel_y, 
                  rel_z) += val;
            }

          }
        }
      }
    }

    //defined in cu. N.B. this _cannot_ be virtual because we are passing the
    //gridmaker to a GPU kernel by value and the vtable is not copied
    template<typename Dtype>
    void setAtomsGPU(unsigned natoms, float4 *ainfos, short *gridindex, 
        qt Q, unsigned ngrids, Dtype *grids);

    template <typename Dtype>
    void setAtomGradientsCPU(const vector<float4>& ainfo, 
        const vector<short>& gridindex, quaternion Q, Dtype* data, 
        vector<float3>& agrad, unsigned offset, unsigned ntypes) {
      unsigned grids_per_dim = (this->dimension - subgrid_dim) / 
        (subgrid_dim + resolution) + 1;
      unsigned ngrids = grids_per_dim * grids_per_dim * grids_per_dim;
      unsigned subgrid_dim_in_points = dim / grids_per_dim;
      boost::multi_array_ref<Dtype, 6> grids(data, boost::extents[ngrids]
          [batch_size][ntypes][subgrid_dim_in_points][subgrid_dim_in_points]
          [subgrid_dim_in_points]);
      zeroAtomGradientsCPU(agrad);
      for (unsigned i = 0, n = ainfo.size(); i < n; ++i) {
        int whichgrid = gridindex[i]; // this is which atom-type channel of the grid to look at
        if (whichgrid >= 0) {
          setAtomGradientCPU(ainfo[i], whichgrid, Q, grids, agrad[i]);
        }
      }
    }

    template<typename Grids>
    void setAtomGradientsCPU(const vector<float4>& ainfo, 
        const vector<short>& gridindex, const quaternion& Q, 
        const Grids& grids, vector<float3>& agrad) {
      zeroAtomGradientsCPU(agrad);
      for (unsigned i = 0, n = ainfo.size(); i < n; ++i) {
        int whichgrid = gridindex[i]; // this is which atom-type channel of the grid to look at
        if (whichgrid >= 0) {
          setAtomGradientCPU(ainfo[i], whichgrid, Q, grids, agrad[i]);
        }
      }
    }

    //TODO: possible to merge this with base version?
    template<typename Grids>
    void setAtomGradientCPU(const float4& ainfo, int whichgrid, 
        const quaternion& Q, const Grids& grids, float3& agrad, 
        bool isrelevance = false) {
      float3 coords;
      if (Q.real() != 0) {//apply rotation
        quaternion p(0, ainfo.x - center.x, ainfo.y - center.y,
            ainfo.z - center.z);
        p = Q * p * (conj(Q) / norm(Q));

        coords.x = p.R_component_2() + center.x;
        coords.y = p.R_component_3() + center.y;
        coords.z = p.R_component_4() + center.z;
      } else {
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
      for (unsigned i = ranges[0].first, iend = ranges[0].second; i < iend;
          ++i) {
        for (unsigned j = ranges[1].first, jend = ranges[1].second;
            j < jend; ++j) {
          for (unsigned k = ranges[2].first, kend = ranges[2].second;
              k < kend; ++k) {
            //convert grid point coordinates to angstroms
            float x = dims[0].x + i * resolution;
            float y = dims[1].x + j * resolution;
            float z = dims[2].x + k * resolution;
            //if we're iterating over subgrids, the LSTM layer requires that the
            //subgrids from different examples in the same batch be interleaved.
            //we need to convert the absolute grid index to the index in this
            //subgrid grid, which is T x B x ntypes x subgrid_dim x subgrid_dim x subgrid_dim
            unsigned grids_per_dim = (dimension - subgrid_dim) / 
              (subgrid_dim + resolution) + 1;
            unsigned subgrid_idx_x = i / (dim / grids_per_dim); 
            unsigned subgrid_idx_y = j / (dim / grids_per_dim); 
            unsigned subgrid_idx_z = k / (dim / grids_per_dim); 
            unsigned rel_x = i % (dim / grids_per_dim); 
            unsigned rel_y = j % (dim / grids_per_dim); 
            unsigned rel_z = k % (dim / grids_per_dim); 
            unsigned grid_idx = (((subgrid_idx_x * grids_per_dim) + subgrid_idx_y) * 
                grids_per_dim + subgrid_idx_z);

            if (isrelevance) 
              accumulateAtomRelevance(coords, radius, x, y, z,
                  grids[grid_idx][batch_idx][whichgrid][rel_x][rel_y][rel_z], 
                  agrad);
            else //true gradient, distance matters
              accumulateAtomGradient(coords, radius, x, y, z,
                  grids[grid_idx][batch_idx][whichgrid][rel_x][rel_y][rel_z], 
                  agrad, whichgrid);
          }
        }
      }
    }

    template<typename Dtype>
    __device__
    void setAtomGradientsGPU(const float4* ainfo, short* gridindices, 
        float3* agrads, const qt Q, const Dtype* grids, 
        unsigned remainder_offset, bool isrelevance = false) {
#ifdef __CUDA_ARCH__
      //TODO: implement
      assert(isrelevance == false);

      int idx = blockDim.x * blockIdx.x + threadIdx.x + remainder_offset;
      int whichgrid = gridindices[idx];
      float4 atom = ainfo[idx];
      float3 coords; 

    if (Q.real() != 0) {//apply rotation
      qt p(atom.x - center.x, atom.y - center.y,
          atom.z - center.z, 0);
      p = Q * p * (Q.conj() / Q.norm());

      coords.x = p.R_component_1() + center.x;
      coords.y = p.R_component_2() + center.y;
      coords.z = p.R_component_3() + center.z;
    } else {
      coords.x = atom.x;
      coords.y = atom.y;
      coords.z = atom.z;
    }

    //get grid index ranges that could possibly be overlapped by atom
    float radius = atom.w;
    float r = radius * radiusmultiple;
    uint2 ranges[3];
    ranges[0] = getrange_gpu(dims[0], coords.x, r);
    ranges[1] = getrange_gpu(dims[1], coords.y, r);
    ranges[2] = getrange_gpu(dims[2], coords.z, r);

    for (unsigned i = ranges[0].x, iend = ranges[0].y; i < iend;
        ++i) {
      for (unsigned j = ranges[1].x, jend = ranges[1].y;
          j < jend; ++j) {
        for (unsigned k = ranges[2].x, kend = ranges[2].y;
            k < kend; ++k) {
          //convert grid point coordinates to angstroms
          float x = dims[0].x + i * resolution;
          float y = dims[1].x + j * resolution;
          float z = dims[2].x + k * resolution;
          unsigned grids_per_dim = (dimension - subgrid_dim) / 
            (subgrid_dim + resolution) + 1;
          unsigned subgrid_idx_x = i / (dim / grids_per_dim); 
          unsigned subgrid_idx_y = j / (dim / grids_per_dim); 
          unsigned subgrid_idx_z = k / (dim / grids_per_dim); 
          unsigned rel_x = i % (dim / grids_per_dim); 
          unsigned rel_y = j % (dim / grids_per_dim); 
          unsigned rel_z = k % (dim / grids_per_dim); 
          unsigned grid_idx = (((subgrid_idx_x * grids_per_dim) + 
                subgrid_idx_y) * grids_per_dim + subgrid_idx_z);

          accumulateAtomGradient(coords, radius, x, y, z, 
                    grids[((((grid_idx * batch_size + batch_idx) * ntypes + 
                          whichgrid) * grids_per_dim + rel_x) * grids_per_dim + 
                      rel_y) * grids_per_dim + rel_z], agrads[idx], whichgrid);
        }
      }
    }
#endif
  }
};

  template <typename Dtype, class GridMakerT>
  __global__ 
  void setAtomGradientGPU(GridMakerT gmaker, 
          float4* ainfo, short* gridindices, float3* agrads, float3 centroid, 
          qt Q, float3 translation, Dtype *diff, int offset, 
          unsigned remainder_offset) {
      gmaker.setAtomGradientsGPU(ainfo, gridindices, agrads, Q, diff + offset, 
          remainder_offset);
  }

template __device__
void GridMaker::setAtomGradientsGPU<double>(const float4* ainfo,
    short* gridindices,
    float3* agrads, const qt Q, const double* grids, unsigned remainder_offset,
    bool isrelevance);
template __device__
void GridMaker::setAtomGradientsGPU<float>(const float4* ainfo,
    short* gridindices,
    float3* agrads, const qt Q, const float* grids, unsigned remainder_offset,
    bool isrelevance);

template __device__
void RNNGridMaker::setAtomGradientsGPU<double>(const float4* ainfo, 
    short* gridindices, float3* agrads, const qt Q, const double* grids, 
    unsigned remainder_offset, bool isrelevance);
template __device__
void RNNGridMaker::setAtomGradientsGPU<float>(const float4* ainfo, 
    short* gridindices, float3* agrads, const qt Q, const float* grids, 
    unsigned remainder_offset, bool isrelevance);

#endif /* _GRIDMAKER_H_ */
