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

#ifndef _POINT_H
#define _POINT_H

#include "basics.h"

namespace T_MESH
{

//! Orientation predicates using filtering on doubles
extern "C" double orient2d(double *, double *, double *);
extern "C" double orient3d(double *, double *, double *, double *);

//! Orientation predicates on PM_Rationals

// orient2D: >0 =0 <0 if (p,q,r) are CCW, aligned, CW respectively
PM_Rational orient2D(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry);


//! Geometric point definition

//! This class represents a point in the Euclidean 3D space. It can be used
//! to  represent  3D vectors originating at (0,0,0) and terminating at the
//! corresponding point. Several methods of  this  class  are  intended  to
//! manipulate  vectors  rather  than  points;  for  example, a call of the
//! method normalize is an actual normalization if the object is a  vector,
//! but  it  has  to  be intended as a projection on the unit sphere if the
//! object is intended to be a point. An object of type Point is a  triplet
//! (x,y,z)  of  coordinates  endowed with a pointer 'info' to possible additional
//! information. Each coordinate is a number of type 'coord' which, by
//! default,  is  a standard double. Operations on points include addition,
//! subtraction, cross and dot product, and many others. This class  implements
//! several useful operations using vector arithmethic. For example,
//! the simple piece of code "A = B*C;" assignes to A the value of the  dot
//! product of B and C.
//! Nearly zero or nearly flat angles are automatically snapped to
//! exactly zero and exactly flat angles if the difference is smaller
//! than the global variable _acos_tolerance. This is the very basic application
//! of our version of the epsilon geometry for robust computation.


class Point
{
 public :
 coord x,y,z;					//!< Coordinates
 void *info;					//!< Further information

 //! Creates a new point with coordinates (0,0,0).
 Point() {x = y = z = 0; info = NULL;}

 //! Creates a new point with the same coordinates as 's'. The info field is not copied.
 Point(const Point *s) {x = s->x; y = s->y; z = s->z; info = NULL;}

 //! Creates a new point with the same coordinates as 's'. The info field is not copied.
 Point(const Point& s) {x = s.x; y = s.y; z = s.z; info = NULL;}

 //! Creates a new point with coordinates (a,b,c).
 Point(const coord& a, const coord& b, const coord& c) {x = a; y = b; z = c; info = NULL;}

 //! Do not remove this. It makes the compiler produce a vtable for this object.
 TMESH_VIRTUAL bool isPoint() const { return true; }

 //! Set the coordinates to (a,b,c).
 void	setValue(const coord& a, const coord& b, const coord& c) {x = a; y = b; z = c;}

 //! Set the coordinates as those of 'p'
 void	setValue(const Point& p) {x = p.x; y = p.y; z = p.z;}

 //! Set the coordinates as those of '*p'
 void	setValue(const Point *p) {x = p->x; y = p->y; z = p->z;}

 //! Returns the vector difference
 Point 	operator-(const Point& p) const {return Point(x-p.x, y-p.y, z-p.z);}

 //! Returns the vector sum
 Point 	operator+(const Point& p) const {return Point(x+p.x, y+p.y, z+p.z);}

 //! Sums another point
 void 	operator+=(const Point& p) {x+=p.x; y+=p.y; z+=p.z;}

 //! Subtracts another point
 void 	operator-=(const Point& p) {x-=p.x; y-=p.y; z-=p.z;}

 //! Returns the Cross Product
 Point 	operator&(const Point& p) const {return Point(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);}

 //! Returns the Dot Product
 coord operator*(const Point& p) const {return (x*p.x+y*p.y+z*p.z);}

 //! Returns the product with a scalar
 Point  operator*(const coord& d) const {return Point(x*d,y*d,z*d);}

 //! Multiplies by a scalar
 void 	operator*=(const coord& m) { x *= m; y *= m; z *= m; }

 //! Divides by a scalar
 void 	operator/=(const coord& m) { x /= m; y /= m; z /= m; }

 //! Returns the vector divided by the scalar
 Point 	operator/(const coord& d) const { return Point(x / d, y / d, z / d); }

 //! TRUE iff coordinates are equal
 bool  	operator==(const Point& p) const {return (x==p.x && y==p.y && z==p.z);}

 //! FALSE iff coordinates are equal
 bool  	operator!=(const Point& p) const {return (x!=p.x || y!=p.y || z!=p.z);}

 //! TRUE iff this is lexycographically smaller than s
 bool operator<(const Point& s) const;

 //! Returns the i'th coordinate
 inline coord& at(unsigned char i) { return (i == 0) ? (x) : ((i == 1) ? (y) : (z)); }

 //! Returns the i'th coordinate
 inline coord& operator[](unsigned char i) { return (i == 0) ? (x) : ((i == 1) ? (y) : (z)); }

 //! Returns the inverse vector
 Point 	inverse() const {return Point(-x,-y,-z);}

 //! Inverts the vector
 void 	invert() {x=-x; y=-y; z=-z;}

 //! TRUE if vector is (0,0,0)
 bool  	isNull() const {return (x==0 && y==0 && z==0);}

 //! Squared distance from origin
 coord squaredLength() const {return (x*x + y*y + z*z);}

 //! Squared distance from '*b'
 coord squaredDistance(const Point *b) const { return (((*(this)) - (*b)).squaredLength()); }

 //! Returns the solution of the linear system Ax = d, where A is a 3x3 matrix whose rows are row1, row2 and row3, d = this
 Point  linearSystem(const Point& row1, const Point& row2, const Point& row3);


 //! Projects the vector on the plane with normal 'n' passing through the origin.
 void   project(const Point *n);

 //! Returns the projection of the point on the straight line though 'a' and 'b'.
 Point  projection(const Point *a, const Point *b) const;

 //! Prints the coordinates of the point to a file handler. stdout is the default.
 void 	printPoint(FILE *fp = stdout) const { fprintf(fp, "%f %f %f,\n", TMESH_TO_FLOAT(x), TMESH_TO_FLOAT(y), TMESH_TO_FLOAT(z)); }		// Debug


 //! Exact orientation test.
 //! Return value is positive iff the tetrahedron (this,a,b,c) has a positive volume;
 //! It is negative iff the tetrahedron (this,a,b,c) has a negative volume;
 //! It is zero iff the tetrahedron (this,a,b,c) has a zero volume.
 coord exactOrientation(const Point *a, const Point *b, const Point *c) const;
 coord	side3D(const Point *p1, const Point *p2, const Point *p3) const { return exactOrientation(p1, p2, p3); }


 //! Exact misalignment test. Returns TRUE iff points are not aligned.
 bool exactMisalignment(const Point *a, const Point *b) const;
 bool 	notAligned(const Point *a, const Point *b) const { return exactMisalignment(a,b); }

 //! Exact planar side test. Returns TRUE iff 'this', Q, A and B are coplanar
 //! and 'this' and Q are (properly) on the same side of A-B.
 //! Warning! Coplanarity is not checked, result is undetermined if
 //! 'this', Q, A and B are not coplanar.
 bool exactSameSideOnPlane(const Point *Q, const Point *A, const Point *B) const;

 //! Itersection point between lines p-q and r-s. Return INFINITE_POINT if lineas are either parallel or degenerate.
 static Point lineLineIntersection(const Point& p, const Point& q, const Point& r, const Point& s);

 //! Itersection point between line p-q and plane r-s-t. Return INFINITE_POINT for parallel/degenerate args.
 static Point linePlaneIntersection(const Point& p, const Point& q, const Point& r, const Point& s, const Point& t);

 //! Squared area of the triangle p-q-r.
 static coord squaredTriangleArea3D(const Point& p, const Point& q, const Point& r);

 //! true if 'p' is a point of the segment v1-v2 (endpoints excluded)
 static bool pointInInnerSegment(const Point *p, const Point *v1, const Point *v2);

 //! true if 'p' is a point of the segment v1-v2 (endpoints included)
 static bool pointInSegment(const Point *p, const Point *v1, const Point *v2);

 //! true if the coplanar point 'p' is in the inner area of v1-v2-v3.
 //! Undetermined if points are not coplanar.
 static bool pointInInnerTriangle(const Point *p, const Point *v1, const Point *v2, const Point *v3);
 
 //! true if the coplanar point 'p' is either in the inner area of v1-v2-v3 or on its border.
 //! Undetermined if points are not coplanar.
 static bool pointInTriangle(const Point *p, const Point *v1, const Point *v2, const Point *v3);

 //! true if (p1-p2) properly intersects (sp1-sp2) at any point (endpoints included).
 //! Collinear overlapping segments are not considered to be properly intersecting.
 static bool segmentsIntersect(const Point *p1, const Point *p2, const Point *sp1, const Point *sp2);

 //! true if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
 //! Collinear overlapping segments are not considered to be properly intersecting.
 static bool innerSegmentsCross(const Point& p1, const Point& p2, const Point& sp1, const Point& sp2);

 //! true if segment (s1-s2) intersects the triangle v1-v2-v3 (border included).
 static bool segmentIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3);

 //! true if segment (s1-s2) intersects the triangle v1-v2-v3 (border included).
 //! Accelerated version - relative orientations are passed as parameters.
 static bool segmentIntersectsTriangle(const Point *s1, const Point *s2, const Point *v1, const Point *v2, const Point *v3, const coord& o1, const coord& o2);




 // FUNCTIONS BELOW THIS LINE MAY RETURN APPROXIMATED/NOT ROBUST RESULTS EVEN WHEN USING RATIONALS



 //! Distance from origin
 double length() const { return sqrt(TMESH_TO_DOUBLE((x*x + y*y + z*z))); }

 //! Divides the vector by its length. If isNull() the application exits with an error.
 void 	normalize();

 //! Rotates the vector around 'axis' by 'ang' radians ccw.
 void  	rotate(const Point& axis, const double& ang);

 //! Distance from 'b'
 double distance(const Point& b) const {return (((*(this))-(b)).length());}

 //! Distance from '*b'
 double distance(const Point *b) const { return (((*(this)) - (*b)).length()); }

 //! Distance from straight line through 'a' and 'b'
 double distanceFromLine(const Point *a, const Point *b) const;

 //! Distance from straight line through 'a' and 'b'. *cc is set to the closest line point.
 double distanceFromLine(const Point *a, const Point *b, Point *cc) const;

 double distanceFromEdge(const Point *a, const Point *b) const; //!< Distance from segment a-b

 //! Distance from segment a-b. *cc is set to the closest edge point.
 double distanceFromEdge(const Point *a, const Point *b, Point *cc) const;

 //! Distance between the straight lines through (this) - l1_p2 and l2_p1 - l2_p2.
 double distanceLineLine(const Point *l1_p2, const Point *l2_p1, const Point *l2_p2) const;

 //!< Angle between this vector and 'v' in radians.
 double getAngle(const Point& v) const;

 //! Angle defined by <a, *this, b> in radians.
 double getAngle(const Point& a, const Point& b) const { return (a - (*this)).getAngle(b - (*this)); }

 //! Angle defined by <*a, *this, *b> in radians.
 double getAngle(const Point *a, const Point *b) const { return ((*a) - (*this)).getAngle((*b) - (*this)); }


 //! Line-line closest point computation.
 //! I SUSPECT THIS CAN BE MADE EXACT...

 //! Computes the closest points of the line passing through this and this2,
 //! and the line passing through p1 and p2. The computed points are used to
 //! initialize the  coordinates  of  cpOnThis  and  cpOnOther.  The  method
 //! returns 0 if the lines are parallel, 1 otherwise.
 int    closestPoints(const Point *this2, const Point *p1, const Point *p2, Point *cpOnThis, Point *cpOnOther) const;
};

//! Lexycographic comparison to be used with jqsort() or abstractHeap.
int xyzCompare(const void *p1, const void *p2);

//! Static point with DBL_MAX coordinates.
extern const Point INFINITE_POINT;

//! Checks whether a point is INFINITE_POINT.
#define IS_FINITE_POINT(p) ((p).x < DBL_MAX)

} //namespace T_MESH

#endif // _POINT_H

