/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2013: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is dual-licensed as follows:                                 *
*                                                                           *
* (1) You may use TMesh as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use TMesh as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#ifndef MARCHING_INTS_H
#define MARCHING_INTS_H

#include "tmesh.h"

namespace T_MESH
{

///////////////////////////////////////////////////////////////
//
// A ray-mesh intersection.
//
///////////////////////////////////////////////////////////////

class mc_ints
{
 public:

 coord ic;		// Distance from ray source
 unsigned char sg;		// 1 if ray enters mesh, 0 if it exits
 ExtVertex *v;	// Used to support polygonization
 Triangle *source; // The triangle that generates this intersection

 mc_ints(coord a, unsigned char b, Triangle *s) { ic = a; sg = b; v = NULL; source = s; }
 ~mc_ints() { if (v) delete(v); }

 static int compare(const void *e1, const void *e2);
};


///////////////////////////////////////////////////////////////
//
// A cubical cell to be polygonized
//
///////////////////////////////////////////////////////////////

class mc_cell
{
 public:
 int x,y,z;			// Coordinates (i.e. cell's position)
 mc_ints *ints[12];	// Intersection at each voxel edge

 mc_cell(int a, int b, int c) {x=a; y=b; z=c;}
 mc_cell(int a, int b, int c, mc_ints *m, int i) : x(a), y(b), z(c)
 {
  for (int n=0; n<12; n++) ints[n] = (n==i)?(m):(NULL);
 }

 void polygonize(Basic_TMesh *tin);
 static int compare(const void *e1, const void *e2);

 void merge(mc_cell *m);

private:
	unsigned char lookup() const;	// The polygonization to be used
	unsigned char lookdown();	// The other polygonization to be used
	void purgeIntersections(); // TRUE if intersections are consistent in the cell
};


///////////////////////////////////////////////////////////////
//
// The marching intersections grid
//
///////////////////////////////////////////////////////////////

class mc_grid
{
 Point origin;		 // Origin for normalization
 coord norm;		 // Normalization factor
protected:													//! < AMF_CHANGE - since T_MESH 2.4-2 >
 Basic_TMesh *tin;		 // Triangulation to remesh				//! < AMF_CHANGE - since T_MESH 2.4-2 >
 List *xy, *xz, *zy; // Axis-parallel rays					//! < AMF_CHANGE - since T_MESH 2.4-2 >
 int numrays;		 // Number of rays per axis				//! < AMF_CHANGE - since T_MESH 2.4-2 >

public:	
 mc_grid(Basic_TMesh *_tin, int n);
 ~mc_grid() {delete [] xy; delete [] xz; delete [] zy;}

 TMESH_VIRTUAL mc_ints * newMcInts(coord a, unsigned char b, Triangle *s){	return new mc_ints(a,b,s);	} //! < AMF_CHANGE - since T_MESH 2.4-2 >

 void remesh(bool simplify_result =false);	
 void simplify();							

 static bool segmentIntersectsTriangle(Point& ev1, Point& ev2, Triangle *t, Point& op);

protected:													//! < AMF_CHANGE - since T_MESH 2.4-2 >
	TMESH_VIRTUAL void sample_triangle(Triangle *t);					//! < AMF_CHANGE - since T_MESH 2.4-2 >
	TMESH_VIRTUAL void createVertices(List *l, int i, int j, int k);	//! < AMF_CHANGE - since T_MESH 2.4-2 >
private:
 void sort();
 void purge();
 void createVertices();
 List *createCells();
 void trackOuterHull();

 static inline coord oceil(const coord& d) { coord c; return ((c = ceil(d)) == d) ? (c + 1) : (c); }
 void purgeList(List *l);
};

} //namespace T_MESH

#endif // MARCHING_INTS_H
