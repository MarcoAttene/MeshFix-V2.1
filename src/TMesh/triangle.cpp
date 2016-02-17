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

#include "triangle.h"
#include <stdlib.h>

namespace T_MESH
{

extern "C" double orient2d(double *, double *, double *);

//////////////////// Constructor //////////////////////

//!< AMF_ADD 1.1>
Triangle::Triangle(){
 mask = 0;
 info = NULL;
}

Triangle::Triangle(Edge *a, Edge *b, Edge *c)
{
 e1 = a;
 e2 = b;
 e3 = c;
 mask = 0;
 info = NULL;
}


//////////////////// Normal vector //////////////////////

Point Triangle::getNormal() const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 Point vd = (((*va)-(*vb))&((*vb)-(*vc)));
 coord l = vd.length();

 if (l == 0) return Point(0,0,0);

 return vd/l;
}


////// Directional vector ////////

Point Triangle::getVector() const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 return (((*va) - (*vb))&((*vb) - (*vc)));
}


/////////////////// Normal consistence check ////////////////////

bool Triangle::checkAdjNor(const Triangle *t) const
{
 Edge *e = commonEdge(t);
 if (e == NULL) return 1;

 Edge *ea = nextEdge(e);
 Edge *eb = t->nextEdge(e);
 if (ea->commonVertex(eb) == ea->commonVertex(e)) return 0;

 return 1;
}


//////////////////////// Triangle area /////////////////////////

double Triangle::area() const
{
 double a = e1->length(), b = e2->length(), c = e3->length();
 if (a==0.0 || b==0.0 || c==0.0) return 0.0;
 double p = (a+b+c)/2.0;
 p = p*(p-a)*(p-b)*(p-c); if (p<0) return 0.0;
 return sqrt(p);
}


/////////////////////// Triangle perimeter /////////////////////

double Triangle::perimeter() const
{
 return e1->length()+e2->length()+e3->length();
}


///////////// Barycenter ///////////////////////

Point Triangle::getCenter() const
{
 Point va = *v1(), vb = *v2(), vc = *v3();
 return (va+vb+vc)/3.0; 
}


///////////////////////// Circlecenter /////////////////////////

Point Triangle::getCircleCenter() const
{
 Point va = *v1(), vb = *v2(), vc = *v3();
 Point q1 = vb-va;
 Point q2 = vc-va;
 Point n = q2&q1;
 Point m1 = e2->getMidPoint();
 Point m2 = e1->getMidPoint();

 return Point(n*va,q1*m1,q2*m2).linearSystem(n,q1,q2);
}


/////// Check wether the point is inside the triangle's bounding ball /////

bool Triangle::inSphere(const Point *p) const
{
 Point c = getCircleCenter();
 coord rad = c.squaredDistance(e1->v1);

 return (p->squaredDistance(&c) < rad);
}


//////////////////// Angle at a vertex /////////////////////

double Triangle::getAngle(const Vertex *v) const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 if (v == va) return v->getAngle(vb, vc);
 if (v == vb) return v->getAngle(va, vc);
 if (v == vc) return v->getAngle(vb, va);

 return -1.0;
}


/////////// Angle between the two directional vectors /////////

double Triangle::getDAngle(const Triangle *t) const
{
 Point thisNormal = getVector();
 Point otherNormal = t->getVector();

 if (thisNormal.isNull() || otherNormal.isNull()) return -1.0;

 return thisNormal.getAngle(otherNormal);
}


///////////// Distance from the plane of the triangle //////////////

double Triangle::distanceFromPoint(const Point *p) const
{
	return sqrt(TMESH_TO_DOUBLE(squaredDistanceFromPoint(p)));
}

///////////// Squared distance from the plane of the triangle //////////////

coord Triangle::squaredDistanceFromPoint(const Point *p) const 
{
 Point CA = e1->toVector()&e2->toVector();
 coord CA2 = CA*CA;

 if (CA2 == 0) return -1.0; 
 coord d = ((CA*(*p))-(CA*(*(e1->v1))));

 return (d*d)/CA2;
}


///////////// Distance of point from the triangle //////////////

double Triangle::pointTriangleDistance(const Point *p, Point *cp) const
{
	return sqrt(TMESH_TO_DOUBLE(pointTriangleSquaredDistance(p)));
}


///////////// Distance of point from the triangle //////////////

coord Triangle::pointTriangleSquaredDistance(const Point *p, Edge **closest_edge, Vertex **closest_vertex) const
{
 Vertex *va = v1(), *vb = v2(), *vc = v3();
 Point n(((*va)-(*vb))&((*vb)-(*vc)));
 if (n.x == 0 && n.y == 0 && n.z == 0) return -1.0;

 coord d1 = ((((*va)-(*vb))&((*vb)-(*p)))*n);
 coord d2 = ((((*vb)-(*vc))&((*vc)-(*p)))*n);
 coord d3 = ((((*vc)-(*va))&((*va)-(*p)))*n);

 if (d1 > 0 && d2 > 0 && d3 > 0) // Closest point in inner triangle
 {
	 if (closest_edge != NULL) *closest_edge = NULL;
	 if (closest_vertex != NULL) *closest_vertex = NULL;
	 return squaredDistanceFromPoint(p);
 }

 if (d2 < 0) { va = vb; vb = vc; if (closest_edge != NULL) *closest_edge = e3; }
 else if (d3 < 0) { vb = va; va = vc; if (closest_edge != NULL) *closest_edge = e1; }
 else if (closest_edge != NULL) *closest_edge = e2;

 Point i(p->projection(va,vb));
 Point p1(i-(*va)); Point p2(i-(*vb));

 if (p1*p2 < 0) // Closest point on interior of one edge
 {
	 return i.squaredDistance(p);
 }

 d1=p1.squaredLength(); d2=p2.squaredLength();
 if (d1 < d2) { if (closest_vertex != NULL) *closest_vertex = va; return p->squaredDistance(va); }
 else { if (closest_vertex != NULL) *closest_vertex = vb; return p->squaredDistance(vb); }
}


/////////// Projection of point 'p' on the plane of the triangle /////

Point Triangle::project(const Point *p) const
{
 Point n = getVector();
 if (n.isNull()) return INFINITE_POINT;
 return Point::linePlaneIntersection(*p, (*p) + n, *(v1()), *(v2()), *(v3()));
}

bool Triangle::isExactlyDegenerate() const
{
	return (!v1()->exactMisalignment((v2()), (v3())));
}

//// get longest edge /////

Edge *Triangle::getLongestEdge() const
{
 coord l1 = e1->squaredLength();
 coord l2 = e2->squaredLength();
 coord l3 = e3->squaredLength();
 if (l1>=l2 && l1>=l3) return e1;
 if (l2>=l1 && l2>=l3) return e2;
 return e3;
}

Edge *Triangle::getCapEdge() const
{
	Edge *e;
	e = e1;  if (Point::pointInInnerSegment(oppositeVertex(e), e->v1, e->v2)) return e;
	e = e2;  if (Point::pointInInnerSegment(oppositeVertex(e), e->v1, e->v2)) return e;
	e = e3;  if (Point::pointInInnerSegment(oppositeVertex(e), e->v1, e->v2)) return e;
	return NULL;
}

///////////// Overlap check ////////////////////

bool Triangle::overlaps() const
{
	return (e1->overlaps() || e2->overlaps() || e3->overlaps());
}

Vertex *Triangle::commonVertex(const Triangle *t2) const
{
	if (hasVertex(t2->v1())) return t2->v1();
	if (hasVertex(t2->v2())) return t2->v2();
	if (hasVertex(t2->v3())) return t2->v3();
	return NULL;
}


/// Debug

void Triangle::printTriangle(FILE *fp) const
{
 v1()->printPoint(fp);
 v2()->printPoint(fp);
 v3()->printPoint(fp);
}


// This can be made more efficient, I guess...

bool Triangle::intersects(const Triangle *t2, bool justproper) const
{
 Vertex *v11, *v12, *v13, *v21, *v22, *v23;

 if (justproper)
 {
	 // This works for non-degenerate triangles. Not sure it will work for degeneracies too.
	 v11 = v1(); v12 = v2(); v13 = v3();
	 v21 = t2->v1(); v22 = t2->v2(); v23 = t2->v3();
		Vertex *eq1 = ((*v11) == (*v21)) ? (v21) : (((*v11) == (*v22)) ? (v22) : (((*v11) == (*v23)) ? (v23) : (NULL)));
		Vertex *eq2 = ((*v12) == (*v21)) ? (v21) : (((*v12) == (*v22)) ? (v22) : (((*v12) == (*v23)) ? (v23) : (NULL)));
		Vertex *eq3 = ((*v13) == (*v21)) ? (v21) : (((*v13) == (*v22)) ? (v22) : (((*v13) == (*v23)) ? (v23) : (NULL)));
		if (eq1 && eq2 && eq3) return false; // Triangles coincide
		Edge *ce1 = NULL, *ce2 = NULL;
		if (eq1 && eq2) { ce1 = e2; ce2 = (t2->e1->hasVertices(eq1, eq2)) ? (t2->e1) : ((t2->e2->hasVertices(eq1, eq2)) ? (t2->e2) : (t2->e3)); }
		if (eq2 && eq3) { ce1 = e3; ce2 = (t2->e1->hasVertices(eq3, eq2)) ? (t2->e1) : ((t2->e2->hasVertices(eq3, eq2)) ? (t2->e2) : (t2->e3)); }
		if (eq3 && eq1) { ce1 = e1; ce2 = (t2->e1->hasVertices(eq3, eq1)) ? (t2->e1) : ((t2->e2->hasVertices(eq3, eq1)) ? (t2->e2) : (t2->e3)); }
		if (ce1)
		{
			Vertex *ov = t2->oppositeVertex(ce2);
			return (ov->exactOrientation(v11, v12, v13) == 0 && ov->exactSameSideOnPlane(oppositeVertex(ce1), ce1->v1, ce1->v2));
		}
		Vertex *cv1 = NULL, *cv2 = NULL;
		if (eq1) { cv1 = v11; cv2 = eq1; }
		if (eq2) { cv1 = v12; cv2 = eq2; }
		if (eq3) { cv1 = v13; cv2 = eq3; }
		if (cv1) // If they share a vertex, intersection occurs if the opposite edge intersect the triangle
		{
			Edge *ee1 = oppositeEdge(cv1), *ee2 = t2->oppositeEdge(cv2);
			return (Point::segmentIntersectsTriangle(ee1->v1, ee1->v2, v21, v22, v23) ||
				Point::segmentIntersectsTriangle(ee2->v1, ee2->v2, v11, v12, v13));
		}
 }
 else
 {
	Edge *ce = commonEdge(t2);
	if (ce) // If they share an edge, intersection occurs only if t1 and t2 overlap
	{
		Vertex *ov = t2->oppositeVertex(ce);
		return (ov->exactOrientation(v1(), v2(), v3()) == 0 && ov->exactSameSideOnPlane(oppositeVertex(ce), ce->v1, ce->v2));
	}

	Vertex *cv = commonVertex(t2);
	v11 = v1(); v12 = v2(); v13 = v3();
	v21 = t2->v1(); v22 = t2->v2(); v23 = t2->v3();
	if (cv) // If they share a vertex, intersection occurs if the opposite edge intersect the triangle
	{
		Edge *ee1 = oppositeEdge(cv), *ee2 = t2->oppositeEdge(cv);
		return (Point::segmentIntersectsTriangle(ee1->v1, ee1->v2, v21, v22, v23) || 
			    Point::segmentIntersectsTriangle(ee2->v1, ee2->v2, v11, v12, v13));
	}
 }

 // Fast reject by bounding box
	coord mx = MIN(v11->x, MIN(v13->x, v12->x));
	if (v21->x < mx && v22->x < mx && v23->x < mx) return false;
	mx = MAX(v11->x, MAX(v13->x, v12->x));
	if (v21->x > mx && v22->x > mx && v23->x > mx) return false;
	mx = MIN(v11->y, MIN(v13->y, v12->y));
	if (v21->y < mx && v22->y < mx && v23->y < mx) return false;
	mx = MAX(v11->y, MAX(v13->y, v12->y));
	if (v21->y > mx && v22->y > mx && v23->y > mx) return false;
	mx = MIN(v11->z, MIN(v13->z, v12->z));
	if (v21->z < mx && v22->z < mx && v23->z < mx) return false;
	mx = MAX(v11->z, MAX(v13->z, v12->z));
	if (v21->z > mx && v22->z > mx && v23->z > mx) return false;

	// Calculate relative orientations
	coord o11 = v11->exactOrientation(v21, v22, v23);
	coord o12 = v12->exactOrientation(v21, v22, v23);
	coord o13 = v13->exactOrientation(v21, v22, v23);
	if ((o11>0 && o12>0 && o13>0) || (o11<0 && o12<0 && o13<0)) return false; // t1 above/below t2
	coord o21 = v21->exactOrientation(v11, v12, v13);
	coord o22 = v22->exactOrientation(v11, v12, v13);
	coord o23 = v23->exactOrientation(v11, v12, v13);
	if ((o21>0 && o22>0 && o23>0) || (o21<0 && o22<0 && o23<0)) return false; // t2 above/below t1

	if (o11 == 0 && o12 == 0 && o13 == 0) // t1 and t2 are coplanar
	{
		if (Point::innerSegmentsCross(v11, v12, v21, v22)) return true;
		if (Point::innerSegmentsCross(v11, v12, v22, v23)) return true;
		if (Point::innerSegmentsCross(v11, v12, v23, v21)) return true;
		if (Point::innerSegmentsCross(v12, v13, v21, v22)) return true;
		if (Point::innerSegmentsCross(v12, v13, v22, v23)) return true;
		if (Point::innerSegmentsCross(v12, v13, v23, v21)) return true;
		if (Point::innerSegmentsCross(v13, v11, v21, v22)) return true;
		if (Point::innerSegmentsCross(v13, v11, v22, v23)) return true;
		if (Point::innerSegmentsCross(v13, v11, v23, v21)) return true;
		return (
			Point::pointInTriangle(v11, v21, v22, v23) ||
			Point::pointInTriangle(v12, v21, v22, v23) ||
			Point::pointInTriangle(v13, v21, v22, v23) ||
			Point::pointInTriangle(v21, v11, v12, v13) ||
			Point::pointInTriangle(v22, v11, v12, v13) ||
			Point::pointInTriangle(v23, v11, v12, v13));
	}
	else return (
		Point::segmentIntersectsTriangle(v11, v12, v21, v22, v23, o11, o12) ||
		Point::segmentIntersectsTriangle(v12, v13, v21, v22, v23, o12, o13) ||
		Point::segmentIntersectsTriangle(v13, v11, v21, v22, v23, o13, o11) ||
		Point::segmentIntersectsTriangle(v21, v22, v11, v12, v13, o21, o22) ||
		Point::segmentIntersectsTriangle(v22, v23, v11, v12, v13, o22, o23) ||
		Point::segmentIntersectsTriangle(v23, v21, v11, v12, v13, o23, o21));
}

} //namespace T_MESH
