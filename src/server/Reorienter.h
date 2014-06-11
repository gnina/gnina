/*
 * ReOrienter.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 *
 *  A class for transforming cartesian coordinates
 */

#ifndef REORIENTER_H_
#define REORIENTER_H_


#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include "common.h"
using namespace Eigen;
using namespace std;


class Reorienter
{
	Matrix3d rotation;
	Vector3d translation;

public:
	Reorienter(): rotation(Matrix3d::Identity()),translation(Vector3d::Zero())
	{

	}

	//modify points by rot/trans
	void reorient(vecv& pts) const
	{
		Vector3d pt;
		for(unsigned i = 0, n = pts.size(); i < n; i++)
		{
			Vector3d pt(pts[i][0],pts[i][1],pts[i][2]);
			pt =  rotation*pt+translation;
			pts[i] = vec(pt[0],pt[1],pt[2]);
		}
	}

	void read(istream& in)
	{
		for(unsigned i = 0; i < 3; i++)
		{
			for(unsigned j = 0; j < 3; j++)
			{
				in >> rotation(i,j);
			}
		}
		for(unsigned i = 0; i < 3; i++)
			in >> translation(i);
	}
};




#endif /* REORIENTER_H_ */
