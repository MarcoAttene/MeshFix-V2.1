/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2012: IMATI-GE / CNR                                         *
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

#include "point.h"
#include <stdlib.h>
#include <limits.h>
#include <errno.h>

namespace T_MESH
{

const Point INFINITE_POINT(DBL_MAX, DBL_MAX, DBL_MAX);


//////// Lexicographic Point comparison //////////

// This can be used with std::sort()
bool Point::operator<(const Point& s) const
{
 if (x<s.x) return true; else if (x>s.x) return false;
 if (y<s.y) return true; else if (y>s.y) return false;
 if (z<s.z) return true; else return false;
}

// This can be used with jqsort
int xyzCompare(const void *a, const void *b)
{
 coord c;

 if ((c=(((Point *)a)->x - ((Point *)b)->x)) < 0) return -1;
 if (c > 0) return 1;
 if ((c=(((Point *)a)->y - ((Point *)b)->y)) < 0) return -1;
 if (c > 0) return 1;
 if ((c=(((Point *)a)->z - ((Point *)b)->z)) < 0) return -1;
 if (c > 0) return 1;

 return 0;
}

//////////////// Normalization /////////////////////////

void Point::normalize()
{
 coord l = length();

 if (l == 0) TMesh::error("normalize : Trying to normalize a null vector !\n");

 x/=l;
 y/=l;
 z/=l;
}


//////////////////// Point rotation ////////////////////
/////////// 'ang' radians CCW around 'axis' ////////////

void Point::rotate(const Point& a, const double& ang)
{
 double l, q[4], m[3][3];
 if ((l = a.length())==0.0) return;
 l = sin(ang/2.0)/l;

 q[0] = TMESH_TO_DOUBLE(a.x)*l;
 q[1] = TMESH_TO_DOUBLE(a.y)*l;
 q[2] = TMESH_TO_DOUBLE(a.z)*l;
 q[3] = cos(ang/2.0);

 m[0][0] = 1.0 - (q[1]*q[1] + q[2]*q[2])*2.0;
 m[0][1] = (q[0] * q[1] + q[2] * q[3])*2.0;
 m[0][2] = (q[2] * q[0] - q[1] * q[3])*2.0;

 m[1][0] = (q[0] * q[1] - q[2] * q[3])*2.0;
 m[1][1] = 1.0 - (q[2] * q[2] + q[0] * q[0])*2.0;
 m[1][2] = (q[1] * q[2] + q[0] * q[3])*2.0;

 m[2][0] = (q[2] * q[0] + q[1] * q[3])*2.0;
 m[2][1] = (q[1] * q[2] - q[0] * q[3])*2.0;
 m[2][2] = 1.0 - (q[1] * q[1] + q[0] * q[0])*2.0;

 q[0] = TMESH_TO_DOUBLE(x); q[1] = TMESH_TO_DOUBLE(y); q[2] = TMESH_TO_DOUBLE(z);
 x = m[0][0]*q[0] + m[1][0]*q[1] + m[2][0]*q[2];
 y = m[0][1]*q[0] + m[1][1]*q[1] + m[2][1]*q[2];
 z = m[0][2]*q[0] + m[1][2]*q[1] + m[2][2]*q[2];
}


///// Project the point on the plane whose normal is 'nor' /////

void Point::project(const Point *nor)
{
 Point pr = (*this)-((*nor)*((*this)*(*nor)));
 x = pr.x; y = pr.y; z = pr.z;
}


////////////// Alignment check /////////////

bool Point::exactMisalignment(const Point *A, const Point *B) const
{
#ifdef USE_HYBRID_KERNEL
	if (coord::use_rationals)
	{
		if (coord::orient2D(x, y, A->x, A->y, B->x, B->y) != 0) return true;
		if (coord::orient2D(y, z, A->y, A->z, B->y, B->z) != 0) return true;
		if (coord::orient2D(z, x, A->z, A->x, B->z, B->x) != 0) return true;
	}
	else
	{
#endif
		double dc[10];

		dc[0] = TMESH_TO_DOUBLE(x); dc[1] = TMESH_TO_DOUBLE(y);
		dc[3] = TMESH_TO_DOUBLE(A->x); dc[4] = TMESH_TO_DOUBLE(A->y);
		dc[6] = TMESH_TO_DOUBLE(B->x); dc[7] = TMESH_TO_DOUBLE(B->y);
		if (orient2d(dc, dc + 3, dc + 6) != 0) return true;

		dc[2] = TMESH_TO_DOUBLE(z); dc[5] = TMESH_TO_DOUBLE(A->z); dc[8] = TMESH_TO_DOUBLE(B->z);
		if (orient2d(dc + 1, dc + 4, dc + 7) != 0) return true;

		dc[9] = dc[6]; dc[6] = dc[3];  dc[3] = dc[0];
		if (orient2d(dc + 2, dc + 5, dc + 8) != 0) return true;
#ifdef USE_HYBRID_KERNEL
	}
#endif

	return false;
}


/////////// Distance from the line passing through A and B ////////

double Point::distanceFromLine(const Point *A, const Point *B) const
{
 Point BA = (*B)-(*A);
 double lba = BA.length();

 if (lba == 0.0) TMesh::error("distanceFromLine : Degenerate line passed !\n");

 return ((((*this)-(*A))&BA).length())/(lba);
}


/////////////////// Distance from a line ///////////////////////
//// 'cc' is initialized as the point of the line whose     ////
//// distance from 'this' is minimum.                       ////

double Point::distanceFromLine(const Point *A, const Point *B, Point *cc) const
{
 Point AB = (*A)-(*B);
 Point AP = (*A)-(*this);
 Point BP = (*B)-(*this);

 if (AP.isNull())
 {
  cc->x = A->x; cc->y = A->y; cc->z = A->z; 
  return 0.0;
 }
 else if (BP.isNull())
 {
  cc->x = B->x; cc->y = B->y; cc->z = B->z; 
  return 0.0;
 }

 coord t = (AB*AB);
 if (t == 0.0) TMesh::error("distanceFromLine : Degenerate line passed !\n");
 else t = (AP*AB)/(-t);
 cc->x = t*AB.x + A->x;
 cc->y = t*AB.y + A->y;
 cc->z = t*AB.z + A->z;
 return distanceFromLine(A,B);
}


////////////// Projection on the line passing through A and B ///////////

Point Point::projection(const Point *A, const Point *B) const
{
 Point BA = (*B)-(*A);
 coord l = BA*BA;
 if (l == 0.0) TMesh::error("projection : Degenerate line passed !\n");

 return ((*A)+(BA*((BA*((*this)-(*A)))/(l))));
}


////////////// Distance from a segment /////////////////

double Point::distanceFromEdge(const Point *A, const Point *B) const
{
 Point AP = (*A)-(*this); double apl = AP.length();
 Point BP = (*B) - (*this); double bpl = BP.length();

 if (apl == 0 || bpl == 0.0) return 0.0;

 Point AB = (*A) - (*B); double abl = AP.length();
 Point BA = (*B)-(*A);

 if (abl*apl == 0.0 || abl*bpl == 0.0) return apl;

 if (AB.getAngle(AP) > PI2) return apl;
 else if (BA.getAngle(BP) > PI2) return bpl;

 return distanceFromLine(A,B);
}

/////////////////// Distance from a segment ///////////////////////
//// 'cc' is initialized as the point of the segment whose     ////
//// distance from 'this' is minimum.                          ////

double Point::distanceFromEdge(const Point *A, const Point *B, Point *cc) const
{
 Point AP = (*A)-(*this); double apl = AP.length();
 Point BP = (*B) - (*this); double bpl = BP.length();

 if (apl == 0) {cc->setValue(A); return 0.0;}
 if (bpl == 0) {cc->setValue(B); return 0.0;}

 Point AB = (*A)-(*B); coord abl = AP.length();
 Point BA = (*B)-(*A);

 if (abl*apl == 0.0 || abl*bpl == 0.0) {cc->setValue(A); return apl;}

 if (AB.getAngle(AP) > PI2) {cc->setValue(A); return apl;}
 else if (BA.getAngle(BP) > PI2) {cc->setValue(B); return bpl;}
 
 coord t = (AB*AB);
 if (t == 0.0) {cc->setValue(A); return apl;}
 else t = (AP*AB)/(-t);
 cc->x = t*AB.x + A->x;
 cc->y = t*AB.y + A->y;
 cc->z = t*AB.z + A->z;
 return distanceFromLine(A,B);
}

///////////////// Angle between two vectors ///////////////

double Point::getAngle(const Point& p) const
{
	return atan2(((*this)&p).length(), TMESH_TO_DOUBLE(((*this)*p)));
}


/////////// Distance of two straight lines ///////////////

double Point::distanceLineLine(const Point *A, const Point *A1, const Point *B1) const
{
 Point uu1 = ((*this)-(*A))&((*A1)-(*B1));
 coord nom = ((*A)-(*A1))*(uu1);
 return FABS(TMESH_TO_DOUBLE(nom)) / (uu1.length());
}


/////////// Solution of a linear system 3 x 3    //////////
///// System Ax = d, where A = (a,b,c) rows, d = this /////

Point Point::linearSystem(const Point& a, const Point& b, const Point& c)
{
 Point ret;
 coord det_A = TMESH_DETERMINANT3X3(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
 if (det_A == 0.0) return INFINITE_POINT;
 ret.x = TMESH_DETERMINANT3X3(x, a.y, a.z, y, b.y, b.z, z, c.y, c.z);
 ret.y = TMESH_DETERMINANT3X3(a.x, x, a.z, b.x, y, b.z, c.x, z, c.z);
 ret.z = TMESH_DETERMINANT3X3(a.x, a.y, x, b.x, b.y, y, c.x, c.y, z);

 return (ret/det_A);
}



//// Computes the closest points of the two lines 'this'-v1 and p1-p2  ////
//// Returns FALSE if the lines are parallel.                          ////

int Point::closestPoints(const Point *v1, const Point *p1, const Point *p2, Point *ptOnThis, Point *ptOnLine2) const
{
 Point pos1 = *this; Point dir1 = (*v1)-pos1;
 Point pos2 = *p1;   Point dir2 = (*p2)-pos2;
 coord d1l = dir1.length(), d2l = dir2.length();

 if (d1l == 0.0 && d2l == 0.0)
  {ptOnThis->setValue(this); ptOnLine2->setValue(p1); return 1;}
 if (d1l*d2l == 0.0)
 {
  if (d1l <= d2l)
   {ptOnThis->setValue(this); distanceFromLine(p1, p2, ptOnLine2); return 1;}
  if (d2l <= d1l)
   {ptOnLine2->setValue(p1); p1->distanceFromLine(this, v1, ptOnThis); return 1;}
 }

 coord ang = dir1.getAngle(dir2);
 if (ang == 0.0 || ang == M_PI) return 0;

 coord s, t, A, B, C, D, E, F, denom;

 denom = ((dir1*dir2)/(d1l*d2l));
 denom = denom*denom - 1;

 dir1.normalize();
 dir2.normalize();

 A = E = dir1*dir2;
 B = dir1*dir1;
 C = (dir1*pos1) - (dir1*pos2);
 D = dir2*dir2;
 F = (dir2*pos1) - (dir2*pos2);

 s = ( C * D - A * F ) / denom;
 t = ( C * E - B * F ) / denom;
 *ptOnThis  = pos1 + (dir1*s);
 *ptOnLine2 = pos2 + (dir2*t);

// Uncomment the following to compute the distance between segments
// if (s < 0 || s > ((*v1)-(*this)).length() || t < 0 || t > ((*p2)-(*p1)).length())
//	return 0;	       // The points does not belong to the edges

 return 1;
}


bool Point::exactSameSideOnPlane(const Point *Q, const Point *A, const Point *B) const
{
#ifdef USE_HYBRID_KERNEL
	if (coord::use_rationals)
	{
		coord o1, o2;
		int s1, s2;

		o1 = coord::orient2D(x, y, A->x, A->y, B->x, B->y);
		o2 = coord::orient2D(Q->x, Q->y, A->x, A->y, B->x, B->y);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(y, z, A->y, A->z, B->y, B->z);
		o2 = coord::orient2D(Q->y, Q->z, A->y, A->z, B->y, B->z);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(z, x, A->z, A->x, B->z, B->x);
		o2 = coord::orient2D(Q->z, Q->x, A->z, A->x, B->z, B->x);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;
	}
	else
	{
#endif
		double dc[13], o1, o2;
		int s1, s2;

		dc[0] = TMESH_TO_DOUBLE(x); dc[1] = TMESH_TO_DOUBLE(y);
		dc[3] = TMESH_TO_DOUBLE(A->x); dc[4] = TMESH_TO_DOUBLE(A->y);
		dc[6] = TMESH_TO_DOUBLE(B->x); dc[7] = TMESH_TO_DOUBLE(B->y);
		dc[9] = TMESH_TO_DOUBLE(Q->x); dc[10] = TMESH_TO_DOUBLE(Q->y);
		o1 = orient2d(dc, dc+3, dc+6);
		o2 = orient2d(dc+9, dc+3, dc+6);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		dc[2] = TMESH_TO_DOUBLE(z); dc[5] = TMESH_TO_DOUBLE(A->z); dc[8] = TMESH_TO_DOUBLE(B->z); dc[11] = TMESH_TO_DOUBLE(Q->z);
		o1 = orient2d(dc+1, dc+4, dc+7);
		o2 = orient2d(dc+10, dc+4, dc+7);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		dc[12] = dc[9]; dc[9] = dc[6]; dc[6] = dc[3];  dc[3] = dc[0];
		o1 = orient2d(dc+2, dc+5, dc+8);
		o2 = orient2d(dc+11, dc+5, dc+8);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;
#ifdef USE_HYBRID_KERNEL
	}
#endif
	return true;
}

coord Point::exactOrientation(const Point *a, const Point *b, const Point *c) const
{
#ifdef USE_HYBRID_KERNEL
	if (coord::use_rationals)
	{
	 return TMESH_DETERMINANT3X3(x - c->x, y - c->y, z - c->z, a->x - c->x, a->y - c->y, a->z - c->z, b->x - c->x, b->y - c->y, b->z - c->z);
	}
	else
	{
#endif
	 double p1[3], p2[3], p3[3], p4[3];
	 p1[0] = TMESH_TO_DOUBLE(x);    p1[1] = TMESH_TO_DOUBLE(y);    p1[2] = TMESH_TO_DOUBLE(z);
	 p2[0] = TMESH_TO_DOUBLE(a->x); p2[1] = TMESH_TO_DOUBLE(a->y); p2[2] = TMESH_TO_DOUBLE(a->z);
	 p3[0] = TMESH_TO_DOUBLE(b->x); p3[1] = TMESH_TO_DOUBLE(b->y); p3[2] = TMESH_TO_DOUBLE(b->z);
	 p4[0] = TMESH_TO_DOUBLE(c->x); p4[1] = TMESH_TO_DOUBLE(c->y); p4[2] = TMESH_TO_DOUBLE(c->z);
	 return orient3d(p1, p2, p3, p4);
#ifdef USE_HYBRID_KERNEL
	}
#endif
}

// Returns the point of intersection between the two lines defined by (p,q) and (r,s) respectively
// Return INFINITE_POINT is lines are parallel or if p==q or r==s
Point Point::lineLineIntersection(const Point& p, const Point& q, const Point& r, const Point& s)
{
	Point da = q - p;
	Point db = s - r;
	Point dc = r - p;
	Point dab = (da&db);

	if (dc*dab != 0.0) return INFINITE_POINT;

	coord k = (((dc&db)*dab) / (dab*dab));
	return p + (da*k);
}

// Returns the point of intersection between the line for (p,q) and the plane for (r,s,t)
// Returns INFINITE_POINT in case of parallelism
Point Point::linePlaneIntersection(const Point& p, const Point& q, const Point& r, const Point& s, const Point& t)
{
	coord den = TMESH_DETERMINANT3X3(p.x - q.x, p.y - q.y, p.z - q.z, s.x - r.x, s.y - r.y, s.z - r.z, t.x - r.x, t.y - r.y, t.z - r.z);
	if (den == 0) return INFINITE_POINT;
	coord num = TMESH_DETERMINANT3X3(p.x - r.x, p.y - r.y, p.z - r.z, s.x - r.x, s.y - r.y, s.z - r.z, t.x - r.x, t.y - r.y, t.z - r.z);
	coord gamma = num / den;
	return p + ((q - p)*gamma);
}

coord Point::squaredTriangleArea3D(const Point& p, const Point& q, const Point& r)
{
	Point pr = (p - r), qr = (q - r);
	Point n = pr&qr;
	return (n*n) / 4;
}


//////////////////////////////////////////////////////////////////
//
// Basic predicates of type 'pointIn'
//
//////////////////////////////////////////////////////////////////

// Returns true if 'p' is a point of the segment v1-v2 (endpoints excluded)
bool Point::pointInInnerSegment(const Point *p, const Point *v1, const Point *v2)
{
	if (!p->exactMisalignment(v1, v2)) // Segment and point aligned
	{
		if (v1->x < v2->x && v1->x < p->x && p->x < v2->x) return true;
		if (v1->y < v2->y && v1->y < p->y && p->y < v2->y) return true;
		if (v1->z < v2->z && v1->z < p->z && p->z < v2->z) return true;
		if (v1->x > v2->x && v1->x > p->x && p->x > v2->x) return true;
		if (v1->y > v2->y && v1->y > p->y && p->y > v2->y) return true;
		if (v1->z > v2->z && v1->z > p->z && p->z > v2->z) return true;
	}
	return false;
}

// Returns true if 'p' is a point of the segment v1-v2 (endpoints included)
bool Point::pointInSegment(const Point *p, const Point *v1, const Point *v2)
{
	return ((*p) == (*(v1)) || (*p) == (*(v2)) || Point::pointInInnerSegment(p, v1, v2));
}

// Returns true if the coplanar point 'p' is in the inner area of 't'.
// Undetermined if p and t are not coplanar.
bool Point::pointInInnerTriangle(const Point *p, const Point *v1, const Point *v2, const Point *v3)
{
	//if (!p->exactSameSideOnPlane(v1, v2, v3)) return false;
	//if (!p->exactSameSideOnPlane(v2, v3, v1)) return false;
	//if (!p->exactSameSideOnPlane(v3, v1, v2)) return false;
	//return true;

	// Less readable, but slightly more efficient (12 predicates instead of 18)

#ifdef USE_HYBRID_KERNEL
	if (coord::use_rationals)
	{
		coord o1, o2, oo2, oo4, oo6;
		int s1, s2;

		o1 = coord::orient2D(p->x, p->y, v2->x, v2->y, v3->x, v3->y);
		o2 = oo2 = coord::orient2D(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->y, p->z, v2->y, v2->z, v3->y, v3->z);
		o2 = oo4 = coord::orient2D(v1->y, v1->z, v2->y, v2->z, v3->y, v3->z);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->z, p->x, v2->z, v2->x, v3->z, v3->x);
		o2 = oo6 = coord::orient2D(v1->z, v1->x, v2->z, v2->x, v3->z, v3->x);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->x, p->y, v3->x, v3->y, v1->x, v1->y);
		o2 = oo2,
			s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->y, p->z, v3->y, v3->z, v1->y, v1->z);
		o2 = oo4;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->z, p->x, v3->z, v3->x, v1->z, v1->x);
		o2 = oo6;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->x, p->y, v1->x, v1->y, v2->x, v2->y);
		o2 = oo2;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->y, p->z, v1->y, v1->z, v2->y, v2->z);
		o2 = oo4;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = coord::orient2D(p->z, p->x, v1->z, v1->x, v2->z, v2->x);
		o2 = oo6;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;
	}
	else
	{
#endif
		int s1, s2;
		double dc[16], o1, o2, oo2, oo4, oo6;

		dc[0] = TMESH_TO_DOUBLE(p->x); dc[1] = TMESH_TO_DOUBLE(p->y); dc[2] = TMESH_TO_DOUBLE(p->z); dc[3] = TMESH_TO_DOUBLE(p->x);
		dc[4] = TMESH_TO_DOUBLE(v1->x); dc[5] = TMESH_TO_DOUBLE(v1->y); dc[6] = TMESH_TO_DOUBLE(v1->z); dc[7] = TMESH_TO_DOUBLE(v1->x);
		dc[8] = TMESH_TO_DOUBLE(v2->x); dc[9] = TMESH_TO_DOUBLE(v2->y); dc[10] = TMESH_TO_DOUBLE(v2->z); dc[11] = TMESH_TO_DOUBLE(v2->x);
		dc[12] = TMESH_TO_DOUBLE(v3->x); dc[13] = TMESH_TO_DOUBLE(v3->y); dc[14] = TMESH_TO_DOUBLE(v3->z); dc[15] = TMESH_TO_DOUBLE(v3->x);

		o1 = orient2d(dc, dc+8, dc+12);
		o2 = oo2 = orient2d(dc+4, dc+8, dc+12);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc+1, dc+9, dc+13);
		o2 = oo4 = orient2d(dc+5, dc+9, dc+13);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc+2, dc+10, dc+14);
		o2 = oo6 = orient2d(dc+6, dc+10, dc+14);
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc, dc+12, dc+4);
		o2 = oo2,
			s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc+1, dc+13, dc+5);
		o2 = oo4;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc+2, dc+14, dc+6);
		o2 = oo6;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc, dc+4, dc+8);
		o2 = oo2;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc+1, dc+5, dc+9);
		o2 = oo4;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;

		o1 = orient2d(dc+2, dc+6, dc+10);
		o2 = oo6;
		s1 = (o1>0) ? (1) : ((o1<0) ? (-1) : (0)); s2 = (o2>0) ? (1) : ((o2<0) ? (-1) : (0));
		if (s1 != s2) return false;
#ifdef USE_HYBRID_KERNEL
	}
#endif
	return true;
}


// Returns true if the coplanar point 'p' is either in the inner area of
// 't' or on its border. Undetermined if p and t are not coplanar.
bool Point::pointInTriangle(const Point *p, const Point *v1, const Point *v2, const Point *v3)
{
	if (Point::pointInSegment(p, v1, v2)) return true;
	else if (Point::pointInSegment(p, v2, v3)) return true;
	else if (Point::pointInSegment(p, v3, v1)) return true;
	else return Point::pointInInnerTriangle(p, v1, v2, v3);
}


//////////////////////////////////////////////////////////////////
//
// Basic predicates of type 'segmentIntersects'
//
//////////////////////////////////////////////////////////////////

// true if (p1-p2) properly intersects (sp1-sp2) at any point (endpoints included).
// Collinear overlapping segments are not considered to be properly intersecting.
bool Point::segmentsIntersect(const Point *p1, const Point *p2, const Point *sp1, const Point *sp2)
{
	return (p1->exactOrientation(p2, sp1, sp2) == 0 && !p1->exactSameSideOnPlane(p2, sp1, sp2) && !sp1->exactSameSideOnPlane(sp2, p1, p2));
}

// Returns true if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
// Collinear overlapping segments are not considered to be properly intersecting.
bool Point::innerSegmentsCross(const Point& p1, const Point& p2, const Point& sp1, const Point& sp2)
{
	if (p1 == sp1 || p1 == sp2 || p2 == sp1 || p2 == sp2) return false;
	return (p1.exactOrientation(&p2, &sp1, &sp2) == 0 && !p1.exactSameSideOnPlane(&p2, &sp1, &sp2) && !sp1.exactSameSideOnPlane(&sp2, &p1, &p2));
}

bool Point::segmentIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3)
{
	coord o1, o2, o3;

	coord mx = MIN(s1->x, s2->x);
	if (v1->x < mx && v2->x < mx && v3->x < mx) return false;
	mx = MAX(s1->x, s2->x);
	if (v1->x > mx && v2->x > mx && v3->x > mx) return false;
	mx = MIN(s1->y, s2->y);
	if (v1->y < mx && v2->y < mx && v3->y < mx) return false;
	mx = MAX(s1->y, s2->y);
	if (v1->y > mx && v2->y > mx && v3->y > mx) return false;
	mx = MIN(s1->z, s2->z);
	if (v1->z < mx && v2->z < mx && v3->z < mx) return false;
	mx = MAX(s1->z, s2->z);
	if (v1->z > mx && v2->z > mx && v3->z > mx) return false;

	o1 = s1->exactOrientation(v1, v2, v3);
	o2 = s2->exactOrientation(v1, v2, v3);
	if (o1 == 0 && o2 == 0)
	{
		if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2)) return true;
		if (Point::pointInInnerTriangle(s1, v1, v2, v3) && Point::pointInInnerTriangle(s2, v1, v2, v3)) return true;
		return false;
	}

	if ((o1>0 && o2>0) || (o1<0 && o2<0)) return false; // s1 and s2 are both above/below v1,v2,v3
	o1 = s1->exactOrientation(s2, v1, v2);
	o2 = s1->exactOrientation(s2, v2, v3);
	if ((o1>0 && o2<0) || (o1<0 && o2>0)) return false;
	o3 = s1->exactOrientation(s2, v3, v1);
	if ((o1>0 && o3<0) || (o1<0 && o3>0)) return false;
	if ((o2>0 && o3<0) || (o2<0 && o3>0)) return false;
	return true;
}

bool Point::segmentIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3, const coord& oo1, const coord& oo2)
{
	// In this case the fast reject by bounding box appears to be a disadvantage ...
	if (oo1 == 0 && oo2 == 0)
	{
		if (!s1->exactSameSideOnPlane(s2, v1, v2) && !v1->exactSameSideOnPlane(v2, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v2, v3) && !v2->exactSameSideOnPlane(v3, s1, s2)) return true;
		if (!s1->exactSameSideOnPlane(s2, v3, v1) && !v3->exactSameSideOnPlane(v1, s1, s2)) return true;
		if (Point::pointInInnerTriangle(s1, v1, v2, v3) && Point::pointInInnerTriangle(s2, v1, v2, v3)) return true;
		return false;
	}

	if ((oo1>0 && oo2>0) || (oo1<0 && oo2<0)) return false; // s1 and s2 are both above/below v1,v2,v3
	coord o1, o2, o3;
	o1 = s1->exactOrientation(s2, v1, v2);
	o2 = s1->exactOrientation(s2, v2, v3);
	if ((o1>0 && o2<0) || (o1<0 && o2>0)) return false;
	o3 = s1->exactOrientation(s2, v3, v1);
	if ((o1>0 && o3<0) || (o1<0 && o3>0)) return false;
	if ((o2>0 && o3<0) || (o2<0 && o3>0)) return false;
	return true;
}

} //namespace T_MESH
