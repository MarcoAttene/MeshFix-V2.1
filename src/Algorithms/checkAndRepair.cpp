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

#include "tmesh.h"
#include "jqsort.h"
#include <stdlib.h>
#include <string.h>

namespace T_MESH
{

////////////// Checks the triangulation's connectivity //////////////
//
// This method should be used when implementing new algorithms to
// check the consistency of the connectivity graph. Because an
// inconsistent graph is not assumed by all the other methods, such
// a flaw is considered critical and the program should terminate.
// If connectivity is ok, NULL is returned, otherwise a string
// describing the error is returned.
//
/////////////////////////////////////////////////////////////////////

const char *Basic_TMesh::checkConnectivity()
{
 Vertex *v;
 Edge *e,*e2;
 Triangle *t;
 Node *n,*m;
 List *ve;

 FOREACHVERTEX(v, n)
 {
  if (v == NULL) return "checkConnectivity: detected NULL element in V list!";
  if (v->e0 == NULL) return "checkConnectivity: detected NULL e0 pointer for a vertex!";
  if (!v->e0->hasVertex(v)) return "checkConnectivity: detected wrong e0 pointer for a vertex!";
 }

 FOREACHEDGE(e, n)
 {
  if (e == NULL) return "checkConnectivity: detected NULL element in E list!";
  if (e->v1 == NULL || e->v2 == NULL) return "checkConnectivity: detected edge with one or two NULL end-points!";
  if (e->v1 == e->v2) return "checkConnectivity: detected edge with two coincident end-points!";
  if (e->t1 == NULL && e->t2 == NULL) return "checkConnectivity: detected edge with no incident triangles!";
  if (e->t1 != NULL)
  {
   if (!e->t1->hasEdge(e)) return "checkConnectivity: detected wrong t1 triangle at an edge";
   if (e->commonVertex(e->t1->nextEdge(e)) == e->v1)
	return "checkConnectivity: Edge orientation does not match t1 normal";
  }
  if (e->t2 != NULL)
  {
   if (!e->t2->hasEdge(e)) return "checkConnectivity: detected wrong t2 triangle at an edge";
   if (e->commonVertex(e->t2->nextEdge(e)) == e->v2)
	return "checkConnectivity: Edge orientation does not match t2 normal";
  }
 }

 FOREACHTRIANGLE(t, n)
 {
  if (t == NULL) return "checkConnectivity: detected NULL element in T list!";
  if (t->e1 == NULL || t->e2 == NULL || t->e3 == NULL) return "checkConnectivity: detected NULL as a triangle edge!";
  if (t->e1 == t->e2 || t->e1 == t->e3 || t->e2 == t->e3) return "checkConnectivity: detected triangle with two coincident edges!";
  if (t->v1() == NULL || t->v2() == NULL || t->v3() == NULL) return "checkConnectivity: triangle edges do not share vertices!";
  if (t->e1->t1 != t && t->e1->t2 != t) return "checkConnectivity: detected triangle with 1st edge not pointing to the triangle itself!";
  if (t->e2->t1 != t && t->e2->t2 != t) return "checkConnectivity: detected triangle with 2nd edge not pointing to the triangle itself!";
  if (t->e3->t1 != t && t->e3->t2 != t) return "checkConnectivity: detected triangle with 3rd edge not pointing to the triangle itself!";
 }

 FOREACHEDGE(e, n)
 {
  ve = e->v1->VE();
  FOREACHVEEDGE(ve, e2, m)
  {
   if (e2 != e && e2->oppositeVertex(e->v1) == e->v2) return "checkConnectivity: detected duplicate edge!";
  }
  if (ve->containsNode(e) == NULL) return "checkConnectivity: detected non manifold vertex!";
  delete(ve);
  ve = e->v2->VE();
  FOREACHVEEDGE(ve, e2, m)
  {
   if (e2 != e && e2->oppositeVertex(e->v2) == e->v1) return "checkConnectivity: detected duplicate edge!";
  }
  if (ve->containsNode(e) == NULL) return "checkConnectivity: detected non manifold vertex!";
  delete(ve);
 }

 return NULL;
}


////////////// Duplicate non-manifold vertices /////////////////////
//
// If a vertex is topologically non-manifold, this data structure
// does not guarantee its functionality. Therefore, in order to use
// the same triangle mesh, this method allows to duplicate such
// vertices. Notice that the data-structure cannot code non-manifold
// edges.
//
////////////////////////////////////////////////////////////////////

int Basic_TMesh::duplicateNonManifoldVertices()
{
 Vertex *v;
 Edge *e, *f;
 Node *n, *m;
 List *ve;
 int dv = 0;

 FOREACHEDGE(e, n)
 {
  ve = e->v1->VE();
  if (ve->containsNode(e) == NULL)
  {
   v = newVertex(e->v1);		//!  
   v->info = e->v1->info;		//! < AMF_CHANGE 1.1-2 > 
   v->mask = 0;					//!  
   V.appendHead(v);
   
   FOREACHVEEDGE(ve, f, m) f->replaceVertex(e->v1, v);
   v->e0 = e->v1->e0;
   e->v1->e0 = e;
   dv++;
  }
  delete(ve);
 }
 FOREACHEDGE(e, n)
 {
  ve = e->v2->VE();
  if (ve->containsNode(e) == NULL)
  {
   v = newVertex(e->v2);		//!  
   v->info = e->v2->info;		//! < AMF_CHANGE 1.1-2 >
   v->mask = 0;					//!
   V.appendHead(v);			     

   FOREACHVEEDGE(ve, f, m) f->replaceVertex(e->v2, v);
   v->e0 = e->v2->e0;
   e->v2->e0 = e;
   dv++;
  }
  delete(ve);
 }

 if (dv) d_boundaries = d_handles = d_shells = 1;

 return dv;
}

////////////// Checks the triangulation geometry //////////////
//// 													 //////
//// Looks for coincident vertices, degenerate triangles //////
//// and overlapping triangles.							 //////
//// If something is wrong returns the closest vertex.   //////
//// 													 //////
///////////////////////////////////////////////////////////////

Vertex *Basic_TMesh::checkGeometry()
{
 int i;
 Vertex *ret = NULL;
 double ang, minda = 0;
 Triangle *t;
 Edge *e;
 Vertex **varr = (Vertex **)V.toArray();
 Edge **evarr;
 Vertex *v1, *v2;
 Node *n;

 if (varr == NULL) TMesh::warning("checkGeometry: Not enough memory. Can't check for coincident vertices.\n");
 else
 {
  jqsort((void **)varr, V.numels(), xyzCompare);
  for (i=0; i<(V.numels()-1); i++)
  {
   v1 = ((Vertex *)varr[i]);
   v2 = ((Vertex *)varr[i+1]);
   if ((*v1)==(*v2))
   {
    ret = v1;
    TMesh::warning("checkGeometry: detected coincident vertices.\n");
    if (v1->getEdge(v2))
    {
     TMesh::warning("               and there is an edge connecting them!\n");
     free(varr);
     return v1;
    }
   }
  }
  free(varr);
 }

 evarr = (Edge **)E.toArray();
 if (evarr == NULL) TMesh::warning("checkGeometry: Not enough memory. Can't check for coincident edges.\n");
 else
 {
  jqsort((void **)evarr, E.numels(), lexEdgeCompare);
  for (i=0; i<(E.numels()-1); i++)
  {
   if (!lexEdgeCompare(evarr[i], evarr[i+1]))
   {
    ret = ((Edge *)evarr[i])->v1;
    TMesh::warning("checkGeometry: detected coincident edges.\n");
   }
  }
  free(evarr);
 }

 FOREACHTRIANGLE(t, n)
 {
  ang = t->getAngle(t->v1());
  if (ang == 0 || ang == M_PI) {TMesh::warning("checkGeometry: degenerate triangle detected.\n"); return t->v1();}
  ang = t->getAngle(t->v2());
  if (ang == 0 || ang == M_PI) {TMesh::warning("checkGeometry: degenerate triangle detected.\n"); return t->v2();}
  ang = t->getAngle(t->v3());
  if (ang == 0 || ang == M_PI) {TMesh::warning("checkGeometry: degenerate triangle detected.\n"); return t->v3();}
 }

 ang = minda = 0;
 FOREACHEDGE(e, n)
  if (e->t1 != NULL && e->t2 != NULL && (ang = e->t1->getDAngle(e->t2)) == M_PI)
   {TMesh::warning("checkGeometry: overlapping triangles detected.\n"); return e->v1;}
  else minda = MAX(minda,ang);
 TMesh::info("checkGeometry: minimum dihedral angle = %f (%f DEGs)\n", M_PI-minda, ((M_PI-minda)*360)/(2*M_PI));
 return ret;
}


///// Merges possible coincident edges //////////

int Basic_TMesh::mergeCoincidentEdges()
{
	V.sort(&xyzCompare);
	Node *n;
	Vertex *v, *pv = (Vertex *)V.head()->data;

	FOREACHVERTEX(v, n) UNMARK_BIT(v, 5);

	Edge *e;
	FOREACHEDGE(e, n) if (e->isOnBoundary()) { MARK_BIT(e->v1, 5); MARK_BIT(e->v2, 5); }

	FOREACHVERTEX(v, n)
	{
		if ((*v) != (*pv) || !IS_BIT(v,5)) pv = v;
		v->info = pv;
		UNMARK_BIT(v, 5);
	}

	// At this point any vertex points (through 'info') to its unique representative (possibly itself)

	FOREACHVERTEX(v, n) v->e0 = NULL;

	FOREACHEDGE(e, n)
	{
		if (e->v1->info != e->v1) e->v1 = (Vertex *)e->v1->info;
		if (e->v2->info != e->v2) e->v2 = (Vertex *)e->v2->info;
		e->v1->e0 = e->v2->e0 = e;
	}
	int rv = removeVertices();

	// At this point the mesh should no longer have duplicated vertices, but may have duplicated edges
	E.sort(&vtxEdgeCompare);
	Edge *pe = (Edge *)E.head()->data;
	FOREACHEDGE(e, n)
	{
		if (!e->isOnBoundary() || vtxEdgeCompare(e, pe)) pe = e;
		e->info = pe;
	}
	FOREACHEDGE(e, n) if (e->info != e)
	{
		Triangle *t1 = e->getBoundaryTriangle();
		Edge *f = ((Edge *)e->info);
		Triangle *t2 = f->getBoundaryTriangle();
		t1->replaceEdge(e, f);
		((f->t1 == NULL) ? (f->t1) : (f->t2)) = t1;
		e->v1 = e->v2 = NULL;
		f->v1->e0 = f->v2->e0 = f;
	}
	removeUnlinkedElements();

	return 1;
}

///// Fix geometric connectivity //////////

bool Basic_TMesh::fixConnectivity(){ //!< AMF_ADD 1.1>
 bool retval = true;
 int i;

 if ((i = removeVertices())) { retval = false; TMesh::warning("%d isolated vertices have been removed.\n", i); }
 if (cutAndStitch()) { retval=false; TMesh::warning("Some cuts were necessary to cope with non manifold configuration.\n"); }
 if (forceNormalConsistence()) { retval = false; TMesh::warning("Some triangles have been reversed to achieve orientation.\n"); }
 if ((i=duplicateNonManifoldVertices())) { retval=false; TMesh::warning("%d non-manifold vertices have been duplicated.\n",i); }
 if ((i=removeDuplicatedTriangles())) { retval=false; TMesh::warning("%d double-triangles have been removed.\n",i); }

 return retval;
}

///// Recompute triangle connectivity //////////

bool Basic_TMesh::rebuildConnectivity(bool fixconnectivity) //!< AMF_CHANGE 1.1>
{
 if (V.numels() == 0) return false;
 V.sort(&xyzCompare);
 Node *n;
 Vertex *v, *pv=(Vertex *)V.head()->data;

 FOREACHVERTEX(v,n)
 {
  if ((*v)!=(*pv)) pv=v;
  v->info = pv;
 }

 // At this point any vertex points (through 'info') to its unique representative (possibly itself)

 FOREACHVERTEX(v,n) v->e0=NULL;

 Edge *e;
 FOREACHEDGE(e, n)
 {
  if (e->v1->info != e->v1) e->v1 = (Vertex *)e->v1->info;
  if (e->v2->info != e->v2) e->v2 = (Vertex *)e->v2->info;
  e->v1->e0 = e->v2->e0 = e;
 }
 int rv = removeVertices();

 // At this point the mesh should no longer have duplicated vertices, but may have duplicated edges

 Triangle *t;
 ExtVertex **var = new ExtVertex *[V.numels()];
 int i=0;
 FOREACHVERTEX(v, n) { v->e0 = NULL; var[i] = new ExtVertex(v); v->info = (void *)i; i++; }
 int nt = T.numels();
 int *triangles = new int[nt*3];
 i = 0; FOREACHTRIANGLE(t, n)
 {
  triangles[i * 3] = (j_voidint)t->v1()->info;
  triangles[i*3+1] = (j_voidint)t->v2()->info;
  triangles[i*3+2] = (j_voidint)t->v3()->info;
  i++;
 }
 T.freeNodes();
 E.freeNodes();
 int v1,v2,v3;
 for (i = 0; i<nt; i++)
 {
  v1 = triangles[i*3];
  v2 = triangles[i*3+1];
  v3 = triangles[i*3+2];
  if (v1!=v2 && v2!=v3 && v1!=v3) CreateIndexedTriangle(var, v1, v2, v3);
 }

 for (i=0; i<V.numels(); i++) delete(var[i]);
 delete var;
 delete [] triangles;

 if(fixconnectivity)	return fixConnectivity();
 else					return true;
}


//////// Eliminates duplicated triangles (i.e. having the same vertices) /////////

int Basic_TMesh::removeDuplicatedTriangles()
{
 Edge *e;
 Node *n;
 Point p;
 int i=0;

 FOREACHEDGE(e, n)
  if (!e->isOnBoundary() && e->t1->oppositeVertex(e) == e->t2->oppositeVertex(e))
  {
	  unlinkTriangle(e->t2);
   i++;
  }
 removeUnlinkedElements();

 if (i)  d_boundaries = d_handles = d_shells = 1;

 return i;
}


//////// Split an edge at all the vertices belonging to its inner segment /////////

int multiSplitEdge(Basic_TMesh *tin, Edge *e)
{
	List splitVertices;
	List triangles, tounmark;
	MARK_BIT(e, 5);
	if (e->t1 != NULL) { triangles.appendTail(e->t1); MARK_BIT(e->t1, 5); }
	if (e->t2 != NULL) { triangles.appendTail(e->t2); MARK_BIT(e->t2, 5); }

	Triangle *t, *y;
	Vertex *v;

	while ((t = (Triangle *)triangles.popHead()) != NULL)
	{
		tounmark.appendHead(t);
		int num_front_edges = 0;
		if (IS_BIT(t->e1, 5)) num_front_edges++;
		if (IS_BIT(t->e2, 5)) num_front_edges++;
		if (IS_BIT(t->e3, 5)) num_front_edges++;
		if (num_front_edges == 3) continue;
		if (num_front_edges == 1)
		{
			Edge *f = (IS_BIT(t->e1, 5)) ? (t->e1) : ((IS_BIT(t->e2, 5)) ? (t->e2) : (t->e3));
			v = t->oppositeVertex(f);
			if (!v->exactMisalignment(e->v1, e->v2))
			{
				if (!IS_BIT(v, 5) && Point::pointInInnerSegment(v, e->v1, e->v2)) { splitVertices.appendTail(v); MARK_BIT(v, 5); }
				y = t->nextEdge(f)->oppositeTriangle(t); if (y != NULL && !IS_BIT(y, 5)) { triangles.appendTail(y); MARK_BIT(y, 5); }
				y = t->prevEdge(f)->oppositeTriangle(t); if (y != NULL && !IS_BIT(y, 5)) { triangles.appendTail(y); MARK_BIT(y, 5); }
			}
		}
		MARK_BIT(t->e1, 5); MARK_BIT(t->e2, 5); MARK_BIT(t->e3, 5);
	}

	Node *n;
	FOREACHVTTRIANGLE((&tounmark), t, n) { UNMARK_BIT(t, 5); UNMARK_BIT(t->e1, 5); UNMARK_BIT(t->e2, 5); UNMARK_BIT(t->e3, 5); }
	FOREACHVVVERTEX((&splitVertices), v, n) UNMARK_BIT(v, 5);

	while (splitVertices.numels())
	{
		coord ad, mind = DBL_MAX;
		Vertex *gv;
		FOREACHVVVERTEX((&splitVertices), v, n) if ((ad = v->squaredDistance(e->v2)) < mind) { gv = v; mind = ad; }
		splitVertices.removeNode(gv);
		tin->splitEdge(e, gv);
	}

	return 1;
}

//////// Split caps and collapse needles to eliminate degenerate triangles /////////

int Basic_TMesh::removeDegenerateTriangles()
{
	Node *n;
	Triangle *t;
	Edge *e, *e1, *e2, *e3, *e4;
	Vertex *ov1, *ov2, *splitvs[2];
	int nov;

	List edges(E);
	FOREACHVEEDGE((&edges), e, n) MARK_BIT(e, 5);

	while ((e = (Edge *)edges.popHead()) != NULL) // Split caps
	{
		UNMARK_BIT(e, 5);
		nov = 0;
		ov1 = (e->t1 != NULL) ? (e->t1->oppositeVertex(e)) : NULL;
		ov2 = (e->t2 != NULL) ? (e->t2->oppositeVertex(e)) : NULL;
		if (ov1 != NULL && Point::pointInInnerSegment(ov1, e->v1, e->v2)) splitvs[nov++] = ov1;
		if (ov2 != NULL && Point::pointInInnerSegment(ov2, e->v1, e->v2)) splitvs[nov++] = ov2;
		if (nov == 1) splitvs[1] = splitvs[0];
		if (nov > 1 && ov1->squaredDistance(e->v1) > ov2->squaredDistance(e->v1)) { splitvs[0] = ov2; splitvs[1] = ov1; }

		if (nov)
		{
			e1 = (e->t1 != NULL) ? (e->t1->nextEdge(e)) : NULL;
			e2 = (e->t1 != NULL) ? (e->t1->prevEdge(e)) : NULL;
			e3 = (e->t2 != NULL) ? (e->t2->nextEdge(e)) : NULL;
			e4 = (e->t2 != NULL) ? (e->t2->prevEdge(e)) : NULL;
			splitEdge(e, splitvs[1]);
			if (nov > 1) splitEdge(e, splitvs[0]);
			e = e1;  if (e != NULL && !IS_BIT(e, 5)) { edges.appendTail(e); MARK_BIT(e, 5); }
			e = e2;  if (e != NULL && !IS_BIT(e, 5)) { edges.appendTail(e); MARK_BIT(e, 5); }
			e = e3;  if (e != NULL && !IS_BIT(e, 5)) { edges.appendTail(e); MARK_BIT(e, 5); }
			e = e4;  if (e != NULL && !IS_BIT(e, 5)) { edges.appendTail(e); MARK_BIT(e, 5); }
		}
	}

	int nc = 0;	// Num of collapses to remove needles

	// Remove needles
	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) == (*e->v2))) if (e->collapse()) nc++;
	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) == (*e->v2)))
	{
		if (e->t1) unlinkTriangle(e->t1);
		if (e->t2) unlinkTriangle(e->t2);
	}
	removeUnlinkedElements();

	int degn = 0;
	FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) degn++;
	if (degn)
	{
		TMesh::warning("removeDegenerateTriangles() - This should not happen!\n");
		FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) MARK_VISIT(t); else UNMARK_VISIT(t);
	}

	return (nc)*((degn) ? (-1) : (1));
}

//int Basic_TMesh::removeDegenerateTriangles()
//{
//	Node *n;
//	Triangle *t;
//	Edge *e;
//
//	// Split caps
//	E.sort(&edgeCompare);
//	Vertex *ov1, *ov2;
//	int nov;
//	bool done;
//
//	do
//	{
//		done = false;
//		FOREACHEDGE(e, n)
//		{
//			nov = 0;
//			if (e->t1 != NULL && Point::pointInInnerSegment((ov1 = e->t1->oppositeVertex(e)), e->v1, e->v2)) nov++;
//			if (e->t2 != NULL && Point::pointInInnerSegment((ov2 = e->t2->oppositeVertex(e)), e->v1, e->v2)) nov += 2;
//			if (nov == 3 && ov1->squaredDistance(e->v1) < ov2->squaredDistance(e->v1)) { splitEdge(e, ov2); splitEdge(e, ov1); }
//			else if (nov == 3 && ov2->squaredDistance(e->v1) < ov1->squaredDistance(e->v1)) { splitEdge(e, ov1); splitEdge(e, ov2); }
//			else if (nov >= 2) splitEdge(e, ov2);
//			else if (nov == 1) splitEdge(e, ov1);
//			if (nov) done = true;
//		}
//	} while (done);
//
//	//FOREACHEDGE(e, n) multiSplitEdge(this, e);
//
//	int nc = 0;	// Num of collapses to remove needles
//
//	// Remove needles
//	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) == (*e->v2))) if (e->collapse()) nc++;
//
//	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) == (*e->v2)))
//	{
//		if (e->t1) unlinkTriangle(e->t1);
//		if (e->t2) unlinkTriangle(e->t2);
//	}
//	removeUnlinkedElements();
//
//	int degn = 0;
//	FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) degn++;
//	if (degn)
//	{
//		FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) MARK_VISIT(t); else UNMARK_VISIT(t);
//	}
//
//	return (nc)*((degn) ? (-1) : (1));
//}

//int Basic_TMesh::removeDegenerateTriangles()
//{
//	Edge *e;
//	Node *n;
//	int nc = 0;	// Num of collapses to remove needles
//	int ns = 0; // Num of swaps to remove internal caps
//	int nr = 0; // Num of removals to remove boundary caps
//
//	// Remove needles
//	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) == (*e->v2))) if (e->collapse()) nc++;
//
//	// Swap or remove caps
//	Point p, w1, w2;
//	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) != (*e->v2)))
//	{
//		if (e->t1 != NULL)
//		{
//			p.setValue(e->t1->oppositeVertex(e));
//			w1 = p - (*e->v1);
//			w2 = p - (*e->v2);
//			if (!w1.isNull() && !w2.isNull() && (w1*w2 < 0) && e->t1->isExactlyDegenerate())
//			{
//				if (e->t2 != NULL) { if (e->swap()) { ns++; continue; } }
//				else { unlinkTriangle(e->t1); nr++; continue; }
//			}
//		}
//		if (e->t2 != NULL)
//		{
//			p.setValue(e->t2->oppositeVertex(e));
//			w1 = p - (*e->v1);
//			w2 = p - (*e->v2);
//			if (!w1.isNull() && !w2.isNull() && (w1*w2 < 0) && e->t2->isExactlyDegenerate())
//			{
//				if (e->t1 != NULL) { if (e->swap()) { ns++; continue; } }
//				else { unlinkTriangle(e->t2); nr++; continue; }
//			}
//		}
//	}
//
//	// Remove possible newly-introduced needles
//	FOREACHEDGE(e, n) if (e->isLinked() && ((*e->v1) == (*e->v2))) if (e->collapse()) nc++;
//
//	removeUnlinkedElements();
//
//	Triangle *t;
//	int degn = 0;
//	FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) degn++;
//	if (degn)
//	{
//		FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) MARK_VISIT(t); else UNMARK_VISIT(t);
//	}
//
//	return (nc + ns + nr)*((degn) ? (-1) : (1));
//}

//int Basic_TMesh::removeDegenerateTriangles()
//{
// Triangle *t;
// Node *n;
// Edge *e;
// int degn = 0, tcs = 0;
//
// if (TMesh::acos_tolerance > 0.0)
// {
//  int collapses, swaps;
//  List todo, *vt;
//
//  do
//  {
//   collapses = swaps = 0;
//
//   FOREACHTRIANGLE(t, n) if (t->isDegenerate())
//    {MARK_BIT(t, 5); todo.appendHead(t);}
//   else UNMARK_BIT(t, 5);
//
//   while (todo.numels())
//   {
//    t = (Triangle *)todo.popHead();
//	UNMARK_BIT(t, 5);
//    if (t->isLinked())
//    {
//     if ((e = t->isCap()) != NULL)
//     {
//      if (e->isOnBoundary()) {unlinkTriangle(t); collapses++;}
//      else if (e->swap())
//      {
//       if (e->t1->overlaps() || e->t2->overlaps() || e->t1->isCap() || e->t2->isCap()) e->swap(1);
//       else
//       {
//        swaps++;
//		if (!IS_BIT(e->t1, 5)) { MARK_BIT(e->t1, 5); todo.appendTail(e->t1); }
//		if (!IS_BIT(e->t2, 5)) { MARK_BIT(e->t2, 5); todo.appendTail(e->t2); }
//       }
//      }
//     }
//     else if ((e = t->isNeedle()) != NULL)
//     {
//      vt = e->v2->VT();
//      if (e->collapse())
//      {
//       collapses++;
//       FOREACHVTTRIANGLE(vt, t, n)
//	   if (t->isLinked() && !IS_BIT(t, 5)) { MARK_BIT(t, 5); todo.appendTail(t); }
//      }
//      else if (!e->isOnBoundary() && e->oppositeTriangle(t)->oppositeVertex(e)->valence() == 3)
//      {
//       if (e->oppositeTriangle(t)->nextEdge(e)->collapse())
//       {
//        MARK_BIT(t, 5); todo.appendHead(t);
//        t = e->oppositeTriangle(t);
//		if (t && !IS_BIT(t, 5)) { MARK_BIT(t, 5); todo.appendTail(t); }
//        collapses++;
//       }
//      }
//      delete(vt);
//     }
//    }
//   }
//
//   if (collapses) removeUnlinkedElements();
//   tcs += collapses; tcs += swaps;
//  } while (collapses+swaps);
//
//  FOREACHTRIANGLE(t, n) if (t->isDegenerate()) degn++;
//
//  if (degn)
//  {
//   FOREACHTRIANGLE(t, n) if (t->isDegenerate()) MARK_VISIT(t); else UNMARK_VISIT(t);
//  }
// }
// else /// This uses exact arithmetics to deduce that triangles are degenerate
// {
//  List triangles;
//  const int MAX_ATTEMPTS = 10;
//
//  FOREACHTRIANGLE(t, n) t->info=0;
//
//  // BIT5 means that the triangle is in the list
//  FOREACHTRIANGLE(t, n)
//  {
//   if (t->isExactlyDegenerate()) {triangles.appendTail(t); MARK_BIT(t, 5);}
//   else UNMARK_BIT(t, 5);
//  }
//
//  while ((t=(Triangle *)triangles.popHead())!=NULL)
//  {
//   UNMARK_BIT(t, 5);
//   if (t->isLinked())
//   {
//	if (t->e1->isDegenerate()) {t->e1->collapse(); tcs++;}
//    else if (t->e2->isDegenerate()) {t->e2->collapse(); tcs++;}
//    else if (t->e3->isDegenerate()) {t->e3->collapse(); tcs++;}
//    else if ((e=t->getLongestEdge())!=NULL)
//    {
//     if (e->swap())
//	 {
//	  tcs++;
//	  t=e->t1;
//	  if (t->isExactlyDegenerate() && !IS_BIT(t, 5) && ((int)t->info < MAX_ATTEMPTS))
//	   {triangles.appendTail(t); MARK_BIT(t, 5); t->info = (void *)(((int)t->info)+1);}
//	  t=e->t2;
//	  if (t->isExactlyDegenerate() && !IS_BIT(t, 5) && ((int)t->info < MAX_ATTEMPTS))
//	   {triangles.appendTail(t); MARK_BIT(t, 5); t->info = (void *)(((int)t->info)+1);}
//	 }
//    }
//   }
//  }
//
//  removeUnlinkedElements();
//
//  FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) degn++;
//  if (degn)
//  {
//   FOREACHTRIANGLE(t, n) if (t->isExactlyDegenerate()) MARK_VISIT(t); else UNMARK_VISIT(t);
//  }
// }
//
// if (degn) tcs=-tcs;
// if (tcs<0) TMesh::info("removeDegenerateTriangles: %d degeneracies could not be removed and have been selected\n",degn);
// return tcs;
//}


bool Basic_TMesh::strongDegeneracyRemoval(int max_iters)
{
 int n, iter_count = 0;
 bool qstatus = TMesh::quiet;

 TMesh::info("Removing degeneracies...\n");
 while ((++iter_count) <= max_iters && removeDegenerateTriangles()<0)
 {
  for (n=1; n<iter_count; n++) growSelection();
  removeSelectedTriangles();
  removeSmallestComponents();
  TMesh::quiet = true; fillSmallBoundaries(E.numels(), false); TMesh::quiet = qstatus;
  coordBackApproximation();
 }

 if (iter_count > max_iters) return false;
 return true;
}

//// If the mesh is made of more than one connected component ////
//// keep only the biggest one and remove all the others.     ////

int Basic_TMesh::removeSmallestComponents()
{
 Node *n,*m;
 List todo;
 List components;
 List *component, *biggest = NULL;
 Triangle *t, *t1, *t2, *t3;
 int nt = 0, gnt = 0;

 if (T.numels() == 0) return 0;

 FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);

 t = ((Triangle *)T.head()->data);
 n = T.head();
 do
 {
  component = new List;
  components.appendHead(component);
  todo.appendHead(t);
  while (todo.numels())
  {
   t = (Triangle *)todo.head()->data;
   todo.removeCell(todo.head());
   if (!IS_BIT(t, 5))
   {
    t1 = t->t1();
    t2 = t->t2();
    t3 = t->t3();

	if (t1 != NULL && !IS_BIT(t1, 5)) todo.appendHead(t1);
	if (t2 != NULL && !IS_BIT(t2, 5)) todo.appendHead(t2);
	if (t3 != NULL && !IS_BIT(t3, 5)) todo.appendHead(t3);

	MARK_BIT(t, 5);
    component->appendTail(t);
   }
  }
  todo.removeNodes();
  for (; n != NULL; n=n->next()) {t = ((Triangle *)n->data); if (!IS_BIT(t, 5)) break;}
 }
 while (n != NULL);

 int num_comps = components.numels();

 FOREACHNODE(components, n)
  if ((nt = ((List *)n->data)->numels()) > gnt) {gnt=nt; biggest = (List *)n->data;}

 FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);

 nt = 0;
 FOREACHNODE(components, n)
  if (((List *)n->data) != biggest)
   FOREACHVTTRIANGLE(((List *)n->data), t, m)
   {
    if (t->e1->v1 != NULL) t->e1->v1->e0 = NULL;
    if (t->e1->v2 != NULL) t->e1->v2->e0 = NULL;
    if (t->e2->v1 != NULL) t->e2->v1->e0 = NULL;
    if (t->e2->v2 != NULL) t->e2->v2->e0 = NULL;
    if (t->e3->v1 != NULL) t->e3->v1->e0 = NULL;
    if (t->e3->v2 != NULL) t->e3->v2->e0 = NULL;
    t->e1->v1 = t->e1->v2 = t->e2->v1 = t->e2->v2 = t->e3->v1 = t->e3->v2 = NULL;
    t->e1 = t->e2 = t->e3 = NULL;
    nt++;
   }

 FOREACHNODE(components, n) delete((List *)n->data);

 if (nt)
 {
  d_boundaries = d_handles = d_shells = 1;
  removeUnlinkedElements();
  return num_comps-1;
 }

 return 0;
}


//// Remove components whose area is < eps_area

int Basic_TMesh::removeSmallestComponents(double eps_area)
{
	Node *n;
	List todo, component;
	Triangle *t, *s;
	int rem_comps=0;
	double pa;

	if (T.numels() == 0) return 0;

	FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);
	t = ((Triangle *)T.head()->data);
	n = T.head();
	do
	{
		todo.appendTail(t); MARK_BIT(t, 5);
		pa = 0.0;
		while ((t=(Triangle *)todo.popHead())!=NULL)
		{
			s = t->t1(); if (s != NULL && !IS_BIT(s, 5)) { todo.appendTail(s); MARK_BIT(s, 5); }
			s = t->t2(); if (s != NULL && !IS_BIT(s, 5)) { todo.appendTail(s); MARK_BIT(s, 5); }
			s = t->t3(); if (s != NULL && !IS_BIT(s, 5)) { todo.appendTail(s); MARK_BIT(s, 5); }
			component.appendTail(t);
			pa += t->area();
		}
		if (pa < eps_area) { rem_comps++; while ((t = (Triangle *)component.popHead()) != NULL) unlinkTriangle(t); }
		else component.removeNodes();

		for (; n != NULL; n = n->next()) { t = ((Triangle *)n->data); if (!IS_BIT(t, 5)) break; }
	} while (n != NULL);

	FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);

	if (rem_comps)
	{
		d_boundaries = d_handles = d_shells = 1;
		removeUnlinkedElements();
	}

	return rem_comps;
}


//// Traverses the triangulation and inverts normals in order ////
//// to make the adjacences consistent.			      ////
//// returns:
//// 0 = mesh was oriented, nothing done
//// 1 = some triangles were flipped to achieve orientation
//// >1 = mesh was not orientable. Cuts were necessary

int Basic_TMesh::forceNormalConsistence()
{
 int ret = 0;
 Node *n;
 Triangle *t;

 FOREACHTRIANGLE(t, n) if (!IS_BIT(t, 5))
  ret |= forceNormalConsistence(t);
 FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);
 return ret;
}

int Basic_TMesh::forceNormalConsistence(Triangle *t0)
{
 Node *n;
 Edge *e;
 List todo, elist;
 Triangle *t, *t1, *t2, *t3;
 int tmp1, tmp2, r=0, wrn = 0, isclosed = 1;

 todo.appendHead(t0);

 while (todo.numels())
 {
  t = (Triangle *)todo.head()->data;
  todo.removeCell(todo.head());
  if (!IS_BIT(t, 5))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
   if (!IS_BIT(t->e1, 5)) { MARK_BIT(t->e1, 5); elist.appendHead(t->e1); }
   if (!IS_BIT(t->e2, 5)) { MARK_BIT(t->e2, 5); elist.appendHead(t->e2); }
   if (!IS_BIT(t->e3, 5)) { MARK_BIT(t->e3, 5); elist.appendHead(t->e3); }

   if (t1 != NULL && !IS_BIT(t1, 5)) { todo.appendHead(t1); if (!t->checkAdjNor(t1)) { t1->invert(); r = 1; } }
   if (t2 != NULL && !IS_BIT(t2, 5)) { todo.appendHead(t2); if (!t->checkAdjNor(t2)) { t2->invert(); r = 1; } }
   if (t3 != NULL && !IS_BIT(t3, 5)) { todo.appendHead(t3); if (!t->checkAdjNor(t3)) { t3->invert(); r = 1; } }

   MARK_BIT(t, 5);
  }
 }

 FOREACHVEEDGE((&(elist)), e, n)
 {
  UNMARK_BIT(e, 5);
  if (isclosed && e->isOnBoundary()) isclosed = 0;
  tmp1 = (e->t1 != NULL)?((e->commonVertex(e->t1->nextEdge(e)) == e->v1)?(-1):(1)):(0);
  tmp2 = (e->t2 != NULL)?((e->commonVertex(e->t2->nextEdge(e)) == e->v2)?(-1):(1)):(0);

  if (tmp1*tmp2 < 0)
  {
   wrn++;
   if (tmp1 == -1) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));
   Edge *ne = newEdge(e->v2, e->v1);
   E.appendHead(ne);
   e->t2->replaceEdge(e, ne);
   ne->t2 = e->t2; e->t2 = NULL;
  } else if (tmp1 == -1 || tmp2 == -1) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));
 }

 if (wrn)
 {
  d_boundaries = d_handles = d_shells = 1;
  TMesh::warning("forceNormalConsistence: Basic_TMesh was not orientable. Cut performed.\n");
 }

 // Though useful in some easy cases, the flip below destroys the orientation
 // when it is set on purpose (e.g. for internal cavities)

 //if (isclosed)
 //{
 // t = topTriangle(t0);
 // if (t->getNormal().z < 0) {flipNormals(t0); r=1;}
 //}

 if (wrn) r |= 2;

 return r;
}


// If possible, swap edges to remove overlaps. When it is not
// enough, remove the overlapping triangles from the mesh.
// return the number of triangles that was necessary to remove.

int Basic_TMesh::removeOverlappingTriangles()
{
 Node *n;
 Edge *e;
 List oved;

 FOREACHEDGE(e, n) if (e->overlaps()) oved.appendHead(e);
 oved.sort(edgeCompare);

 for (n=oved.tail(); n!=NULL; n=n->prev())
 {
  e = (Edge *)n->data;
  if (e->overlaps() && e->swap())
  {
   if (e->t1->isExactlyDegenerate() || e->t2->isExactlyDegenerate()) {e->swap(1); continue;}
   if (e->t1->nextEdge(e)->overlaps()) { e->swap(1); continue; }
   if (e->t1->prevEdge(e)->overlaps()) { e->swap(1); continue; }
   if (e->t2->nextEdge(e)->overlaps()) { e->swap(1); continue; }
   if (e->t2->prevEdge(e)->overlaps()) { e->swap(1); continue; }
  }
 }

 int nr = 0;
 for (n=oved.tail(); n!=NULL; n=n->prev())
 {
  e = (Edge *)n->data;
  if (e->overlaps()) { unlinkTriangle(e->t1); unlinkTriangle(e->t2); nr++; }
 }
 if (nr)
 {
	 removeUnlinkedElements();
	 d_boundaries = d_handles = d_shells = 1;
 }
 
 return nr*2;
}


bool Basic_TMesh::meshclean(int max_iters, int inner_loops)
{
 bool ni, nd;
 Triangle *t;
 Node *m;

 deselectTriangles();
 invertSelection();

 for (int n=0; n<max_iters; n++)
 {
  TMesh::info("********* ITERATION %d *********\n",n);
  nd = strongDegeneracyRemoval(inner_loops);
  deselectTriangles(); invertSelection();
  ni = strongIntersectionRemoval(inner_loops);
  if (ni && nd)
  {
   FOREACHTRIANGLE(t, m) if (t->isExactlyDegenerate()) ni=false;
   if (ni) return true;
  }
 }

 return false;
}

} //namespace T_MESH
