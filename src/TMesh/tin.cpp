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

#include "tin.h"
#include <stdlib.h>
#include <string.h>

namespace T_MESH
{

	Vertex *	Basic_TMesh::newVertex(){								return new Vertex();						}	//!< AMF_ADD 1.1>
	Vertex *	Basic_TMesh::newVertex(const coord &x, const coord &y, const coord &z){ return new Vertex(x, y, z); }	//!< AMF_ADD 1.1>
	Vertex *	Basic_TMesh::newVertex(Point *p){						return new Vertex(p);						}	//!< AMF_ADD 1.1>	
	Vertex *	Basic_TMesh::newVertex(Point &p){						return new Vertex(p);						}	//!< AMF_ADD 1.1>
	Vertex *	Basic_TMesh::newVertex(Vertex *v){						return new Vertex(v);						}	//!< AMF_ADD 1.1-2>
	Edge *		Basic_TMesh::newEdge(Vertex *s, Vertex *d){				return new Edge(s, d);						}	//!< AMF_ADD 1.1>
	Edge *		Basic_TMesh::newEdge(Edge *e){							return new Edge(e->v1,e->v2);				}	//!< AMF_ADD 1.1-2>
	Triangle *	Basic_TMesh::newTriangle(){								return new Triangle();						}	//!< AMF_ADD 1.1>
	Triangle *	Basic_TMesh::newTriangle(Edge *a, Edge *b, Edge *c){	return new Triangle(a, b, c);				}	//!< AMF_ADD 1.1>

//////////////////////////////////////////////////////////////////
//                                                              //
//    C L A S S   C O N S T R U C T O R S                       //
//                                                              //
//////////////////////////////////////////////////////////////////


///////////////////// Constructor (Empty) ////////////////////

Basic_TMesh::Basic_TMesh()
{
 info=NULL;
 n_boundaries = n_handles = n_shells = 0;
 d_boundaries = d_handles = d_shells = 0;
}


//////////////////// Constructor (Pre-defined) ///////////////

Basic_TMesh::Basic_TMesh(const char *tin_definition) { init(tin_definition); }

void Basic_TMesh::init(const char *tin_definition)
{
 info=NULL;
 if (!strcmp(tin_definition, "triangle"))
 {
  Vertex *v1 = newVertex(0,0,0);
  Vertex *v2 = newVertex(2,0,0);
  Vertex *v3 = newVertex(1,1,0);
  Edge   *e1 = newEdge(v1,v2); v1->e0 = e1;
  Edge   *e2 = newEdge(v2,v3); v2->e0 = e2;
  Edge   *e3 = newEdge(v3,v1); v3->e0 = e3;
  Triangle *t1 = newTriangle(e1,e2,e3);
  e1->t1 = t1; e1->t2 = NULL;
  e2->t1 = t1; e2->t2 = NULL;
  e3->t1 = t1; e3->t2 = NULL;
  V.appendHead(v1); V.appendHead(v2); V.appendHead(v3);
  T.appendHead(t1);
  E.appendHead(e1); E.appendHead(e2); E.appendHead(e3);

  n_boundaries = 1;
  n_handles = 0;
  n_shells = 1;
  d_boundaries = d_handles = d_shells = 0;
 }
 else if (!strcmp(tin_definition, "tetrahedron"))
 {
  Vertex *v1 = newVertex(-1,-1.4142136,0);
  Vertex *v2 = newVertex(-1,1.4142136,0);
  Vertex *v3 = newVertex(1,0,-1.4142136);
  Vertex *v4 = newVertex(1,0,1.4142136);
  Edge   *e1 = newEdge(v1,v2); v1->e0 = e1;
  Edge   *e2 = newEdge(v2,v3); v2->e0 = e2;
  Edge   *e3 = newEdge(v3,v1); v3->e0 = e3;
  Edge   *e4 = newEdge(v1,v4); v4->e0 = e4;
  Edge   *e5 = newEdge(v2,v4);
  Edge   *e6 = newEdge(v3,v4);
  Triangle *t1 = newTriangle(e1,e2,e3);
  Triangle *t2 = newTriangle(e1,e4,e5);
  Triangle *t3 = newTriangle(e2,e5,e6);
  Triangle *t4 = newTriangle(e3,e6,e4);
  e1->t1 = t1; e1->t2 = t2;
  e2->t1 = t1; e2->t2 = t3;
  e3->t1 = t1; e3->t2 = t4;
  e4->t1 = t2; e4->t2 = t4;
  e5->t1 = t3; e5->t2 = t2;
  e6->t1 = t4; e6->t2 = t3;
  V.appendHead(v1); V.appendHead(v2); V.appendHead(v3); V.appendHead(v4);
  T.appendHead(t1); T.appendHead(t2); T.appendHead(t3); T.appendHead(t4);
  E.appendHead(e1); E.appendHead(e2); E.appendHead(e3);
  E.appendHead(e4); E.appendHead(e5); E.appendHead(e6);

  n_boundaries = 0;
  n_handles = 0;
  n_shells = 1;
  d_boundaries = d_handles = d_shells = 0;
 }
 else if (!strcmp(tin_definition, "cube"))
 {
  const double crds[8][3] = {{0, 0, 0},{1, 0, 0},{1, 1, 0},{0, 1, 0},{0, 0, 1},{1, 0, 1},{1, 1, 1},{0, 1, 1}};
  const int tris[12][3] = {{3, 2, 1},{3, 1, 0},{4, 5, 6},{4, 6, 7},{7, 6, 2},{7, 2, 3},{0, 1, 5},{0, 5, 4},{1, 2, 6},{1, 6, 5},{3, 0, 4},{3, 4, 7}};
  ExtVertex *ev[8];
  for (int i=0; i<8; i++)
  {
   Vertex *v = newVertex(crds[i][0], crds[i][1], crds[i][2]);
   ev[i] = new ExtVertex(v);
   V.appendTail(v);
  }
  for (int i=0; i<12; i++)
   CreateIndexedTriangle(ev, tris[i][0], tris[i][1], tris[i][2]);
  for (int i=0; i<8; i++) delete ev[i];

  n_boundaries = 1;
  n_handles = 0;
  n_shells = 1;
  d_boundaries = d_handles = d_shells = 0;
 }
 else if (!strcmp(tin_definition, "cylinder"))
 {
  const double crds[8][2] = {{1,0},{0.7,0.7},{0,1},{-0.7,0.7},{-1,0},{-0.7,-0.7},{0,-1},{0.7,-0.7}};
  ExtVertex *ev[8];
  for (int i=0; i<16; i++)
  {
   Vertex *v = newVertex(crds[i%8][0], crds[i%8][1], (i<8)?(-1):(1));
   ev[i] = new ExtVertex(v);
   V.appendTail(v);
  }

  for (int i=0; i<8; i++)
  {
   CreateIndexedTriangle(ev, i, (i+1)%8, i+8);
   CreateIndexedTriangle(ev, i+8, (i+1)%8, 8+(i+1)%8);
  }
  for (int i=0; i<6; i++)
  {
   CreateIndexedTriangle(ev, 0, i+2, i+1);
   CreateIndexedTriangle(ev, 8, i+9, i+10);
  }
  for (int i=0; i<8; i++) delete ev[i];

  n_boundaries = 1;
  n_handles = 0;
  n_shells = 1;
  d_boundaries = d_handles = d_shells = 0;
 }
 else TMesh::error("Unknown triangulation type '%s'\n",tin_definition);
}


///////////////////// Cloning TIN ///////////////////////////

Basic_TMesh::Basic_TMesh(const Basic_TMesh *tin, const bool clone_info) { init(tin, clone_info); }

void Basic_TMesh::init(const Basic_TMesh *tin, const bool clone_info)
{
 info=NULL;
 Node *n;
 Vertex *v, *nv;
 Edge *e, *ne;
 Triangle *t, *nt;

 int i;
 void **t_info = new void *[tin->T.numels()];
 i=0; FOREACHVTTRIANGLE((&(tin->T)), t, n) t_info[i++]=t->info;
 void **e_info = new void *[tin->E.numels()];
 i=0; FOREACHVEEDGE((&(tin->E)), e, n) e_info[i++]=e->info;
 void **v_info = new void *[tin->V.numels()];
 i=0; FOREACHVVVERTEX((&(tin->V)), v, n) v_info[i++]=v->info;

 FOREACHVVVERTEX((&(tin->V)), v, n)
  {nv=newVertex(v); V.appendTail(nv); v->info = nv;}

 FOREACHVEEDGE((&(tin->E)), e, n)
  {ne=newEdge((Vertex *)e->v1->info, (Vertex *)e->v2->info); E.appendTail(ne); e->info = ne;}

 FOREACHVTTRIANGLE((&(tin->T)), t, n)
  {nt=newTriangle((Edge *)t->e1->info,(Edge *)t->e2->info,(Edge *)t->e3->info); T.appendTail(nt); t->info = nt;}

 FOREACHVVVERTEX((&(tin->V)), v, n) {((Vertex *)v->info)->e0 = (Edge *)v->e0->info; v->info = NULL;}

 FOREACHVEEDGE((&(tin->E)), e, n)
  {((Edge *)e->info)->t1 = (e->t1)?((Triangle *)e->t1->info):(NULL); ((Edge *)e->info)->t2 = (e->t2)?((Triangle *)e->t2->info):(NULL); e->info = NULL;}

 i=0; FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info=t_info[i++];
 i=0; FOREACHVEEDGE((&(tin->E)), e, n) e->info=e_info[i++];
 i=0; FOREACHVVVERTEX((&(tin->V)), v, n) v->info=v_info[i++];

 if (clone_info)
 {
  i=0; FOREACHTRIANGLE(t, n) t->info=t_info[i++];
  i=0; FOREACHEDGE(e, n) e->info=e_info[i++];
  i=0; FOREACHVERTEX(v, n) v->info=v_info[i++];  
 }
 delete [] t_info; delete [] e_info; delete [] v_info;

 d_boundaries = d_handles = d_shells = 1;
}


//// Creates a new Basic_TMesh out of a connected component of an existing Basic_TMesh.
//// If 'keep_reference' is set to 'true', each element of the existing mesh keeps a
//// pointer to the corresponding new element in the 'info' field.

Basic_TMesh::Basic_TMesh(const Triangle *t0, const bool keep_reference) { init(t0, keep_reference); }

void Basic_TMesh::init(const Triangle *t0, const bool keep_reference)
{
 info=NULL;
 List todo(t0), st, sv, se;
 Node *n;
 Triangle *t, *nt;
 Vertex *v, *nv;
 Edge *e, *ne;

 t=(Triangle *)t0;
 MARK_VISIT2(t);

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  st.appendHead(t);
  nt=t->t1(); if (nt != NULL && !IS_VISITED2(nt)) {MARK_VISIT2(nt); todo.appendHead(nt);}
  nt=t->t2(); if (nt != NULL && !IS_VISITED2(nt)) {MARK_VISIT2(nt); todo.appendHead(nt);}
  nt=t->t3(); if (nt != NULL && !IS_VISITED2(nt)) {MARK_VISIT2(nt); todo.appendHead(nt);}
 }

 FOREACHVTTRIANGLE((&st), t, n)
 {
  UNMARK_VISIT2(t);
  e = t->e1; if (!IS_VISITED2(e)) {MARK_VISIT2(e); se.appendHead(e);}
  e = t->e2; if (!IS_VISITED2(e)) {MARK_VISIT2(e); se.appendHead(e);}
  e = t->e3; if (!IS_VISITED2(e)) {MARK_VISIT2(e); se.appendHead(e);}
  v = t->v1(); if (!IS_VISITED2(v)) {MARK_VISIT2(v); sv.appendHead(v);}
  v = t->v2(); if (!IS_VISITED2(v)) {MARK_VISIT2(v); sv.appendHead(v);}
  v = t->v3(); if (!IS_VISITED2(v)) {MARK_VISIT2(v); sv.appendHead(v);}
 }

 FOREACHVVVERTEX((&sv), v, n)
  {UNMARK_VISIT2(v); nv=newVertex(v); V.appendTail(nv); v->info = nv;}

 FOREACHVEEDGE((&se), e, n)
  {UNMARK_VISIT2(e); ne=newEdge((Vertex *)e->v1->info, (Vertex *)e->v2->info); E.appendTail(ne); e->info = ne;}

 FOREACHVTTRIANGLE((&st), t, n)
  {nt=newTriangle((Edge *)t->e1->info,(Edge *)t->e2->info,(Edge *)t->e3->info); T.appendTail(nt); t->info = nt;}

 FOREACHVVVERTEX((&sv), v, n) ((Vertex *)v->info)->e0 = (Edge *)v->e0->info;

 FOREACHVEEDGE((&se), e, n)
  {((Edge *)e->info)->t1 = (e->t1)?((Triangle *)e->t1->info):(NULL); ((Edge *)e->info)->t2 = (e->t2)?((Triangle *)e->t2->info):(NULL);}

 if (!keep_reference)
 {
  FOREACHVVVERTEX((&sv), v, n) v->info = NULL;
  FOREACHVEEDGE((&se), e, n) e->info = NULL;
  FOREACHVTTRIANGLE((&st), t, n) t->info = NULL;
 }

 eulerUpdate();
}


Basic_TMesh *Basic_TMesh::split()
{
 deselectTriangles();
 Triangle *t = (Triangle *)T.head()->data;
 selectConnectedComponent(t);
 Basic_TMesh *stin = createSubMeshFromSelection(t);
 removeSelectedTriangles();
 return stin;
}

///////////////////// Destructor ///////////////////////////

Basic_TMesh::~Basic_TMesh()
{
 T.freeNodes();
 V.freeNodes();
 E.freeNodes();
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    P R I M I T I V E   C O N S T R U C T I O N               //
//                                                              //
//////////////////////////////////////////////////////////////////


//////////////////// Creates an edge ////////////////////////////

Edge *Basic_TMesh::CreateEdge(Vertex *v1, Vertex *v2)
{
 Edge *e;
 
 if ((e = v1->getEdge(v2)) != NULL) return e;
 
 e = newEdge(v1,v2);
 v1->e0 = e;
 v2->e0 = e;
 E.appendHead(e);

 return e;
}


//////////////////// Creates an edge ////////////////////////////

Edge *Basic_TMesh::CreateEdge(ExtVertex *v1, ExtVertex *v2, const bool check)
{
 Edge *e;
 Node *n;

 if (check)
  FOREACHVEEDGE((&(v1->VE)), e, n)
   if (e->oppositeVertex(v1->v) == v2->v) return e;

 e = newEdge(v1->v,v2->v);
 if (v1->v->e0 == NULL) v1->v->e0 = e;
 if (v2->v->e0 == NULL) v2->v->e0 = e;
 v1->VE.appendHead(e);
 v2->VE.appendHead(e);
 E.appendHead(e);

 return e;
}


///////////////////// Creates a triangle //////////////////////////

Triangle *Basic_TMesh::CreateTriangle(Edge *e1, Edge *e2, Edge *e3)
{
 Triangle *tt, **at1, **at2, **at3;
 
 if (e1->commonVertex(e2) == e1->v2 && e1->t1 == NULL) at1 = &(e1->t1);
 else if (e1->commonVertex(e2) == e1->v1 && e1->t2 == NULL) at1 = &(e1->t2);
 else return NULL;
 if (e2->commonVertex(e3) == e2->v2 && e2->t1 == NULL) at2 = &(e2->t1);
 else if (e2->commonVertex(e3) == e2->v1 && e2->t2 == NULL) at2 = &(e2->t2);
 else return NULL;
 if (e3->commonVertex(e1) == e3->v2 && e3->t1 == NULL) at3 = &(e3->t1);
 else if (e3->commonVertex(e1) == e3->v1 && e3->t2 == NULL) at3 = &(e3->t2);
 else return NULL;

 tt = newTriangle(e1,e2,e3);
 *at1 = *at2 = *at3 = tt;
 T.appendHead(tt);

 MARK_VISIT(tt);				// < AMF_ADD - since IMATI-STL 2.4-1 >

 d_boundaries = d_handles = d_shells = 1;

 return tt;
}


///////////////////// Creates an unoriented triangle //////////////////////////

Triangle *Basic_TMesh::CreateUnorientedTriangle(Edge *e1, Edge *e2, Edge *e3)
{
 Triangle *tt, **at1, **at2, **at3;
 
 if (e1->t1 == NULL) at1 = &(e1->t1);
 else if (e1->t2 == NULL) at1 = &(e1->t2);
 else return NULL;
 if (e2->t1 == NULL) at2 = &(e2->t1);
 else if (e2->t2 == NULL) at2 = &(e2->t2);
 else return NULL;
 if (e3->t1 == NULL) at3 = &(e3->t1);
 else if (e3->t2 == NULL) at3 = &(e3->t2);
 else return NULL;

 tt = newTriangle(e1,e2,e3);
 *at1 = *at2 = *at3 = tt;
 T.appendHead(tt);

 return tt;
}


////////////// Euler operatior: Create edge and triangle //////////////////////

Triangle *Basic_TMesh::EulerEdgeTriangle(Edge *e2, Edge *e3)
{
 Vertex *cv = e2->commonVertex(e3);
 Triangle *adj = (e2->t1 == NULL)?(e2->t2):(e2->t1);
 if (cv == NULL || !e2->isOnBoundary() || !e3->isOnBoundary()) return NULL;

 Edge *e1 = CreateEdge(e2->oppositeVertex(cv), e3->oppositeVertex(cv));
 if (adj->nextEdge(e2)->hasVertex(cv)) return CreateTriangle(e1,e3,e2);
 return CreateTriangle(e1,e2,e3);
}


/////////// Creation of two triangles bridging two boundary edges /////////////

Edge *Basic_TMesh::bridgeBoundaries(Edge *gve, Edge *gwe)
{
	if (gve == gwe || !gve->isOnBoundary() || !gwe->isOnBoundary()) return NULL;

	Triangle *t;
	Vertex *v = gve->commonVertex(gwe);
	if (v != NULL)
	{
		t = EulerEdgeTriangle(gve, gwe);
		return gve;
	}

	Vertex *gv = (gve->t1) ? (gve->v1) : (gve->v2);
	Vertex *gw = (gwe->t1) ? (gwe->v2) : (gwe->v1);
	Vertex *gvn = gve->oppositeVertex(gv);
	Vertex *gwn = gwe->oppositeVertex(gw);

	Edge *je = CreateEdge(gv, gw);
	Edge *je2 = CreateEdge(gwn, gvn);
	Edge *je1 = CreateEdge(gv, gwn);

	t = CreateTriangle(je, gwe, je1);
	t = CreateTriangle(je1, je2, gve);

	return je1;
}

//////////////////////////////////////////////////////////////////
//                                                              //
//    P R I M I T I V E   D E S T R U C T I O N                 //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Unlinks the triangle (elements are not removed from the lists) //////////

void Basic_TMesh::unlinkTriangle(Triangle *t)
{
 Vertex *v1 = t->v1(), *v2 = t->v2(), *v3 = t->v3();
 Edge *e1 = t->e1, *e2 = t->e2, *e3 = t->e3;

 int v1nm = (v1->isOnBoundary() && !e1->isOnBoundary() && !e2->isOnBoundary());
 int v2nm = (v2->isOnBoundary() && !e2->isOnBoundary() && !e3->isOnBoundary());
 int v3nm = (v3->isOnBoundary() && !e3->isOnBoundary() && !e1->isOnBoundary());

 v1->e0 = ((e2->isOnBoundary())?(e1):(e2));
 v2->e0 = ((e3->isOnBoundary())?(e2):(e3));
 v3->e0 = ((e1->isOnBoundary())?(e3):(e1));

 e1->replaceTriangle(t, NULL);
 e2->replaceTriangle(t, NULL);
 e3->replaceTriangle(t, NULL);

 if (e1->isIsolated() && e2->isIsolated()) v1->e0 = NULL;
 if (e2->isIsolated() && e3->isIsolated()) v2->e0 = NULL;
 if (e3->isIsolated() && e1->isIsolated()) v3->e0 = NULL;
 if (e1->isIsolated()) e1->v1 = e1->v2 = NULL;
 if (e2->isIsolated()) e2->v1 = e2->v2 = NULL;
 if (e3->isIsolated()) e3->v1 = e3->v2 = NULL;
 t->e1 = t->e2 = t->e3 = NULL;

 Vertex *nv;
 Edge *e;
 List *ve;
 Node *n;

 if (v1nm)
 {
  nv = newVertex(v1->x, v1->y, v1->z);
  nv->e0 = v1->e0;
  ve = v1->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(v1, nv);
  delete(ve);
  v1->e0 = e1;
  V.appendHead(nv);
 }
 if (v2nm)
 {
  nv = newVertex(v2->x, v2->y, v2->z);
  nv->e0 = v2->e0;
  ve = v2->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(v2, nv);
  delete(ve);
  v2->e0 = e2;
  V.appendHead(nv);
 }
 if (v3nm)
 {
  nv = newVertex(v3->x, v3->y, v3->z);
  nv->e0 = v3->e0;
  ve = v3->VE();
  FOREACHVEEDGE(ve, e, n) e->replaceVertex(v3, nv);
  delete(ve);
  v3->e0 = e3;
  V.appendHead(nv);
 }
}


//// Unlinks the triangle (elements are not removed from the lists)              ////
//// Differently from the above method, non-manifold vertices are not duplicated ////

void Basic_TMesh::unlinkTriangleNoManifold(Triangle *t)
{
 Edge *e1 = t->e1, *e2 = t->e2, *e3 = t->e3;

 e1->replaceTriangle(t, NULL);
 e2->replaceTriangle(t, NULL);
 e3->replaceTriangle(t, NULL);

 if (e1->isIsolated()) e1->v1 = e1->v2 = NULL;
 if (e2->isIsolated()) e2->v1 = e2->v2 = NULL;
 if (e3->isIsolated()) e3->v1 = e3->v2 = NULL;
 t->e1 = t->e2 = t->e3 = NULL;
}


///// Removes all the triangles with NULL edges /////

int Basic_TMesh::removeTriangles()
{
 Node *n;
 Triangle *t;
 int r = 0;

 n = T.head();
 while (n != NULL)
 {
  t = (Triangle *)n->data;
  n = n->next();
  if (t->e1 == NULL || t->e2 == NULL || t->e3 == NULL)
  {
   r++;
   T.removeCell((n!=NULL)?(n->prev()):T.tail());
   delete t;
  }
 }

 d_boundaries = d_handles = d_shells = 1;

 return r;
}


///// Removes all the edges with NULL vertices /////

int Basic_TMesh::removeEdges()
{
 Node *n;
 Edge *e;
 int r = 0;

 n = E.head();
 while (n != NULL)
 {
  e = (Edge *)n->data;
  n = n->next();
  if (e->v1 == NULL || e->v2 == NULL)
  {
   r++;
   E.removeCell((n!=NULL)?(n->prev()):E.tail());
   delete e;
  }
 }

 d_boundaries = d_handles = d_shells = 1;

 return r;
}


/////////// Removes all the vertices with e0 field = NULL ////////////

int Basic_TMesh::removeVertices()
{
 Node *n;
 Vertex *v;
 int r = 0;

 n = V.head();
 while (n != NULL)
 {
  v = (Vertex *)n->data;
  n = n->next();
  if (v->e0 == NULL)
  {
   r++;
   V.removeCell((n!=NULL)?(n->prev()):V.tail());
   delete v;
  }
 }

 d_boundaries = d_handles = d_shells = 1;

 return r;
}



//////////////////////////////////////////////////////////////////
//                                                              //
//    S E L E C T I O N   M A N A G E M E N T                   //
//                                                              //
//////////////////////////////////////////////////////////////////


//////////////// Deselect all the triangles /////////////

void Basic_TMesh::deselectTriangles()
{
 Triangle *t;
 Node *n;
 FOREACHTRIANGLE(t, n) UNMARK_VISIT(t);
}


//////// Removes all the selected (IS_VISITED) triangles /////////////

void Basic_TMesh::removeSelectedTriangles()
{
 Node *n;
 Triangle *t;

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) unlinkTriangle(t);
 removeUnlinkedElements();
}


// Mark all the triangles having at least one boundary vertex as 'selected'

int Basic_TMesh::selectBoundaryTriangles()
{
 Node *n;
 Edge *e;
 Vertex *v, *v1, *v2, *v3;
 Triangle *t;
 int ns=0;

 FOREACHEDGE(e, n) if (e->isOnBoundary()) {MARK_VISIT(e->v1); MARK_VISIT(e->v2);}

 FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3)) {MARK_VISIT(t); ns++;}
 }
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 return ns;
}


// Grows the current selection (1 triangle width)

int Basic_TMesh::growSelection()
{
 Node *n;
 Vertex *v, *v1, *v2, *v3;
 Triangle *t;
 int ns=0;

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  MARK_VISIT(v1); MARK_VISIT(v2); MARK_VISIT(v3);
 }
 FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3)) {MARK_VISIT(t); ns++;}
 }

 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 return ns;
}


// Shrinks the current selection (1 triangle width)

void Basic_TMesh::shrinkSelection()
{
 Node *n;
 Vertex *v, *v1, *v2, *v3;
 Triangle *t;

 FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  MARK_VISIT(v1); MARK_VISIT(v2); MARK_VISIT(v3);
 }
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3)) UNMARK_VISIT(t);
 }

 FOREACHVERTEX(v, n) UNMARK_VISIT(v);
}


// Toggles the selection status of the triangles

void Basic_TMesh::invertSelection(Triangle *t0)
{
 Node *n;
 Triangle *t;

 if (t0 != NULL)
 {
  List totoggle(t0);
  Triangle *s;
  bool unmark = IS_VISITED(t0);
  if (unmark) UNMARK_VISIT(t0); else MARK_VISIT(t0);
  while ((t = (Triangle *)totoggle.popHead()) != NULL)
  {
   if ((s = t->t1()) != NULL && ((IS_VISITED(s) && unmark) || (!IS_VISITED(s) && !unmark)))
    {if (unmark) UNMARK_VISIT(s); else MARK_VISIT(s); totoggle.appendTail(s);}
   if ((s = t->t2()) != NULL && ((IS_VISITED(s) && unmark) || (!IS_VISITED(s) && !unmark)))
    {if (unmark) UNMARK_VISIT(s); else MARK_VISIT(s); totoggle.appendTail(s);}
   if ((s = t->t3()) != NULL && ((IS_VISITED(s) && unmark) || (!IS_VISITED(s) && !unmark)))
    {if (unmark) UNMARK_VISIT(s); else MARK_VISIT(s); totoggle.appendTail(s);}
  }
 }
 else
 {
  FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) UNMARK_VISIT(t); else MARK_VISIT(t);
 }
}


void Basic_TMesh::reselectSelection(Triangle *t0)
{
 if (!IS_VISITED(t0)) return;

 Node *n;
 Triangle *t, *s;
 List triList(t0);
 MARK_VISIT2(t0);

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  if ((s = t->t1()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t2()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
  if ((s = t->t3()) != NULL && !IS_VISITED2(s) && IS_VISITED(s)) {triList.appendHead(s); MARK_VISIT2(s);}
 }

 FOREACHTRIANGLE(t, n) if (!IS_VISITED2(t)) UNMARK_VISIT(t); else UNMARK_VISIT2(t);
}

// Creates a new mesh out of a selection.

Basic_TMesh *Basic_TMesh::createSubMeshFromSelection(Triangle *t0, bool keep_ref)
{
 Triangle *t,*s, *nt;
 Node *n;

 if (t0 != NULL && !IS_VISITED(t0)) return NULL;

 Basic_TMesh *tin = newObject();
 Vertex *v,*nv;
 Edge *e, *ne;
 List triList, sT, sE, sV;

 if (t0 != NULL)
 {
  triList.appendHead(t0); MARK_BIT(t0,5);
  while (triList.numels())
  {
	  t = (Triangle *)triList.popHead();
	  sT.appendHead(t);
	  if ((s = t->t1()) != NULL && !IS_BIT(s, 5) && IS_VISITED(s)) { triList.appendHead(s); MARK_BIT(s, 5); }
	  if ((s = t->t2()) != NULL && !IS_BIT(s, 5) && IS_VISITED(s)) { triList.appendHead(s); MARK_BIT(s, 5); }
	  if ((s = t->t3()) != NULL && !IS_BIT(s, 5) && IS_VISITED(s)) { triList.appendHead(s); MARK_BIT(s, 5); }
  }
  FOREACHVTTRIANGLE((&sT), t, n)
  {
	  UNMARK_BIT(t->e1, 5); UNMARK_BIT(t->e2, 5); UNMARK_BIT(t->e3, 5);
	  UNMARK_BIT(t->v1(), 5); UNMARK_BIT(t->v2(), 5); UNMARK_BIT(t->v3(), 5);
  }
  FOREACHVTTRIANGLE((&sT), t, n)
  {
	  if (!IS_BIT(t->e1, 5)) { sE.appendHead(t->e1); MARK_BIT(t->e1, 5); }
	  if (!IS_BIT(t->e2, 5)) { sE.appendHead(t->e2); MARK_BIT(t->e2, 5); }
	  if (!IS_BIT(t->e3, 5)) { sE.appendHead(t->e3); MARK_BIT(t->e3, 5); }
	  if ((v = t->v1()) && !IS_BIT(v, 5)) { sV.appendHead(v); MARK_BIT(v, 5); }
	  if ((v = t->v2()) && !IS_BIT(v, 5)) { sV.appendHead(v); MARK_BIT(v, 5); }
	  if ((v = t->v3()) && !IS_BIT(v, 5)) { sV.appendHead(v); MARK_BIT(v, 5); }
  }
 }
 else
 {
  FOREACHEDGE(e, n) UNMARK_BIT(e, 5);
  FOREACHVERTEX(v, n) UNMARK_BIT(v, 5);
  FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
   {sT.appendHead(t); MARK_BIT(t->e1, 5); MARK_BIT(t->e2, 5); MARK_BIT(t->e3, 5);}
  FOREACHEDGE(e, n) if (IS_BIT(e,5))
   {sE.appendHead(e); MARK_BIT(e->v1, 5); MARK_BIT(e->v2, 5);}
  FOREACHVERTEX(v, n) if (IS_BIT(v,5)) sV.appendHead(v);
 }

 FOREACHVEEDGE((&sE), e, n) e->v1->e0 = e->v2->e0 = e;

 int i;
 void **v_info = NULL, **e_info = NULL, **t_info = NULL;
 if (!keep_ref)
 {
  v_info = new void *[sV.numels()];
  i=0; FOREACHVVVERTEX((&sV), v, n) v_info[i++] = v->info;
  e_info = new void *[sE.numels()];
  i=0; FOREACHVEEDGE((&sE), e, n) e_info[i++] = e->info;
  t_info = new void *[sT.numels()];
  i=0; FOREACHVTTRIANGLE((&sT), t, n) t_info[i++] = t->info;
 }

 FOREACHVVVERTEX((&sV), v, n)
  {nv=newVertex(v); tin->V.appendTail(nv); v->info = nv;}

 FOREACHVEEDGE((&sE), e, n)
  {ne=newEdge((Vertex *)e->v1->info, (Vertex *)e->v2->info); tin->E.appendTail(ne); e->info = ne;}

 FOREACHVTTRIANGLE((&sT), t, n)
  {nt=newTriangle((Edge *)t->e1->info,(Edge *)t->e2->info,(Edge *)t->e3->info); tin->T.appendTail(nt); t->info = nt; nt->info = t;}

 FOREACHVVVERTEX((&sV), v, n) ((Vertex *)v->info)->e0 = (Edge *)v->e0->info;

 FOREACHVEEDGE((&sE), e, n)
 {
  ((Edge *)e->info)->t1 = (e->t1 && IS_VISITED(e->t1))?((Triangle *)e->t1->info):(NULL);
  ((Edge *)e->info)->t2 = (e->t2 && IS_VISITED(e->t2))?((Triangle *)e->t2->info):(NULL);
 }

 i=0; if (!keep_ref) FOREACHVVVERTEX((&sV), v, n) v->info = v_info[i++];
 i=0; if (!keep_ref) FOREACHVEEDGE((&sE), e, n) e->info = e_info[i++];
 i=0; if (!keep_ref) FOREACHVTTRIANGLE((&sT), t, n) t->info = t_info[i++];

 FOREACHVTTRIANGLE((&sT), t, n) UNMARK_BIT(t, 5);
 FOREACHVEEDGE((&sE), e, n) UNMARK_BIT(e, 5);
 FOREACHVVVERTEX((&sV), v, n) UNMARK_BIT(v, 5);

 if (!sT.numels()) {delete(tin); return NULL;}

 tin->duplicateNonManifoldVertices();
 tin->eulerUpdate();

 return tin;
}

Basic_TMesh *Basic_TMesh::createSubMeshFromTriangle(Triangle *t0)
{
	Basic_TMesh *ntin = newObject("triangle");

	Vertex *v1 = (Vertex *)ntin->V.head()->data;
	Vertex *v2 = (Vertex *)ntin->V.head()->next()->data;
	Vertex *v3 = (Vertex *)ntin->V.head()->next()->next()->data;
	v1->setValue(t0->v1());
	v2->setValue(t0->v3());
	v3->setValue(t0->v2());
	((Triangle *)ntin->T.head()->data)->info = t0->info;

	return ntin;
}

///// Marks all the triangles within distance L as selected //////

int Basic_TMesh::selectSphericalRegion(Triangle *t, const double L, const Point *center)
{
 List *reg = getRegion(t, L, center);
 Node *n;
 Triangle *s;
 int nt=0;

 FOREACHVTTRIANGLE(reg, s, n) {MARK_VISIT(s); nt++;}
 delete(reg);

 return nt;
}


///// Deselects all the triangles within distance L //////

int Basic_TMesh::deselectSphericalRegion(Triangle *t, const double L, const Point *center)
{
 List *reg = getRegion(t, L, center);
 Node *n;
 Triangle *s;
 int nt=0;

 FOREACHVTTRIANGLE(reg, s, n) {UNMARK_VISIT(s); nt++;}
 delete(reg);

 return nt;
}


///// Selects all the triangles within distance L which were already //////
///// selected. Deselects the others.                                //////

void Basic_TMesh::reselectSphericalRegion(Triangle *t, const double L, const Point *center)
{
 List *reg = getRegion(t, L, center);
 Node *n;
 Triangle *s;

 FOREACHVTTRIANGLE(reg, s, n) MARK_VISIT2(s);
 FOREACHTRIANGLE(s, n) if (IS_VISITED(s) && !IS_VISITED2(s)) UNMARK_VISIT(s);
 FOREACHVTTRIANGLE(reg, s, n) UNMARK_VISIT2(s);
 delete(reg);
}


//// Remove all the selected triangles and re-triangulate /////

bool Basic_TMesh::retriangulateSelectedRegion()
{
 List ttbr;
 Node *n;
 Triangle *u;
 Point nor;
 FOREACHTRIANGLE(u, n) if (IS_VISITED(u))
  {ttbr.appendHead(u); nor = nor+(u->getNormal()*u->area());}

 if (ttbr.numels() < 2)
 {
  TMesh::warning("retriangulateRegion: Nothing to retriangulate.\n");
  return 0;
 }

 FOREACHVTTRIANGLE((&(ttbr)), u, n)
  if (u->getNormal()*nor <= 0.0)
  {
   TMesh::warning("retriangulateRegion: Too complex geometry. Can't retriangulate.\n");
   return 0;
  }

 if (!isSelectionSimple(&ttbr))
 {
  TMesh::warning("retriangulateRegion: Non-simple region. Can't retriangulate.\n");
  return 0;
 }

 List *ms = getRegionInternalVertices(&ttbr);

 FOREACHVTTRIANGLE((&(ttbr)), u, n) unlinkTriangle(u);
 Edge *e = ((Edge *)ms->head()->data);
 List *vl = ((List *)ms->head()->next()->data);
 TriangulateHole(e, vl);
 delete(vl);
 delete(ms);
 removeUnlinkedElements();

 return 1;
}


////////// Check wether 's' represents a simple selection ////////

bool Basic_TMesh::isSelectionSimple(List *s)
{
 if (!s->numels()) return 0; // Empty region is not simple

 Node *n;
 Triangle *ta, *t = (Triangle *)s->head()->data;
 List bdr, top(t);
 MARK_VISIT2(t);
 int nv=0;

 while (top.numels())
 {
  t = (Triangle *)top.popHead();
  nv++;
  ta=t->t1(); if (ta && IS_VISITED(ta) && !IS_VISITED2(ta)) {MARK_VISIT2(ta); top.appendHead(ta);}
  else if (ta == NULL) break; else if (!IS_VISITED(ta)) bdr.appendHead(t->e1);
  ta=t->t2(); if (ta && IS_VISITED(ta) && !IS_VISITED2(ta)) {MARK_VISIT2(ta); top.appendHead(ta);}
  else if (ta == NULL) break; else if (!IS_VISITED(ta)) bdr.appendHead(t->e2);
  ta=t->t3(); if (ta && IS_VISITED(ta) && !IS_VISITED2(ta)) {MARK_VISIT2(ta); top.appendHead(ta);}
  else if (ta == NULL) break; else if (!IS_VISITED(ta)) bdr.appendHead(t->e3);
 }

 FOREACHVTTRIANGLE(s, t, n) UNMARK_VISIT2(t);
 if (top.numels()) return 0; // Mesh-boundary in selection
 if (nv != s->numels()) return 0; // Disconnected selection

 Edge *e, *f, *ge=NULL, *e0;
 List *ve;
 FOREACHVEEDGE((&(bdr)), e, n) MARK_VISIT(e);
 int nae;

 nv = 0;
 e = e0 = (Edge *)bdr.head()->data;
 Vertex *v = e->v1;

 do
 {
  nv++;
  v = e->oppositeVertex(v);
  ve = v->VE();
  nae=0; FOREACHVEEDGE(ve, f, n) if (f!=e && IS_VISITED(f)) {ge=f; nae++;}
  delete(ve);
  if (nae > 1) break;
  e=ge;
 } while (e != e0);

 FOREACHVEEDGE((&(bdr)), e, n) UNMARK_VISIT(e);
 if (nv != bdr.numels()) return 0; // Non-simple selection

 return 1;
}


//// Unmarks all the elements of the triangulation ////
//// but leaves selected triangles marked.         ////

void Basic_TMesh::unmarkEverythingButSelections()
{
 Vertex *v;
 Edge *e;
 Triangle *t;
 Node *n;
 FOREACHVERTEX(v, n) v->mask = 0;
 FOREACHEDGE(e, n) e->mask = 0;
 FOREACHTRIANGLE(t, n) t->mask &= (unsigned char)1;
}


//// Selects all the triangles of the shell containing 't' ////

int Basic_TMesh::selectConnectedComponent(Triangle *t0, bool sos)
{
 List todo;
 Triangle *t, *t1, *t2, *t3;
 int ns = 0;

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (!IS_VISITED(t))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && !IS_VISITED(t1) && (!(sos && IS_SHARPEDGE(t->e1)))) todo.appendHead(t1);
   if (t2 != NULL && !IS_VISITED(t2) && (!(sos && IS_SHARPEDGE(t->e2)))) todo.appendHead(t2);
   if (t3 != NULL && !IS_VISITED(t3) && (!(sos && IS_SHARPEDGE(t->e3)))) todo.appendHead(t3);

   MARK_VISIT(t); ns++;
  }
 }

 return ns;
}


//// Deselects all the triangles of the shell containing 't' ////

int Basic_TMesh::deselectConnectedComponent(Triangle *t0, bool sos)
{
 List todo;
 Triangle *t, *t1, *t2, *t3;
 int ns = 0;

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (IS_VISITED(t))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && IS_VISITED(t1) && (!(sos && IS_SHARPEDGE(t->e1)))) todo.appendHead(t1);
   if (t2 != NULL && IS_VISITED(t2) && (!(sos && IS_SHARPEDGE(t->e2)))) todo.appendHead(t2);
   if (t3 != NULL && IS_VISITED(t3) && (!(sos && IS_SHARPEDGE(t->e3)))) todo.appendHead(t3);

   UNMARK_VISIT(t); ns++;
  }
 }

 return ns;
}


// Append to the current mesh a copy of all the elements of 'src'.
// The newly created elements form a new selection.

void Basic_TMesh::append(Basic_TMesh *src)
{
 deselectTriangles();
 Basic_TMesh cb(src);
 cb.invertSelection();
 V.joinTailList(&(cb.V));
 E.joinTailList(&(cb.E));
 T.joinTailList(&(cb.T));
 d_boundaries = d_handles = d_shells = 1;
}


// Move all the elements of 't' to this mesh and delete 't' itself.
void Basic_TMesh::moveMeshElements(Basic_TMesh *t, bool delInput)
{
	V.joinTailList(&(t->V));
	E.joinTailList(&(t->E));
	T.joinTailList(&(t->T));
	d_boundaries = d_handles = d_shells = 1;
	if(delInput)	delete t;					
}

//////////////////////////////////////////////////////////////////
//                                                              //
//    R E G I O N   M A N I P U L A T I O N                     //
//                                                              //
//////////////////////////////////////////////////////////////////

///// Make a list with all the triangles within distance L //////

List *Basic_TMesh::getRegion(Triangle *t, const double L, const Point *center)
{
 List triList, *toRemove = new List;
 if (t->v1()->distance(center) > L) return toRemove;
 if (t->v2()->distance(center) > L) return toRemove;
 if (t->v3()->distance(center) > L) return toRemove;

 Triangle *s;
 Node *n;

 triList.appendHead(t);
 MARK_BIT(t,3);

 while(triList.numels() > 0)
 {
  t = (Triangle *)triList.head()->data;
  triList.removeCell(triList.head());
  toRemove->appendHead(t);

  if ((s = t->t1()) != NULL && !IS_BIT(s,3) && s->oppositeVertex(t->e1)->distance(center) <= L)
   {triList.appendHead(s); MARK_BIT(s,3);}
  if ((s = t->t2()) != NULL && !IS_BIT(s,3) && s->oppositeVertex(t->e2)->distance(center) <= L)
   {triList.appendHead(s); MARK_BIT(s,3);}
  if ((s = t->t3()) != NULL && !IS_BIT(s,3) && s->oppositeVertex(t->e3)->distance(center) <= L)
   {triList.appendHead(s); MARK_BIT(s,3);}
 }

 FOREACHVTTRIANGLE(toRemove, s, n) UNMARK_BIT(s, 3);

 return toRemove;
}


///// Unlink all the triangles within distance L //////

void Basic_TMesh::removeRegion(Triangle *t, const double L, const Point *center)
{
 List triList, toRemove;
 Node *n;
 Triangle *s;

 triList.appendHead(t);
 MARK_VISIT(t);

 while(triList.numels() > 0)
 {
  t = (Triangle *)triList.head()->data;
  triList.removeCell(triList.head());
  toRemove.appendHead(t);

  if ((s = t->t1()) != NULL && !IS_VISITED(s) && s->oppositeVertex(t->e1)->distance(center) <= L)
   {triList.appendHead(s); MARK_VISIT(s);}
  if ((s = t->t2()) != NULL && !IS_VISITED(s) && s->oppositeVertex(t->e2)->distance(center) <= L)
   {triList.appendHead(s); MARK_VISIT(s);}
  if ((s = t->t3()) != NULL && !IS_VISITED(s) && s->oppositeVertex(t->e3)->distance(center) <= L)
   {triList.appendHead(s); MARK_VISIT(s);}
 }

 for (n = toRemove.tail(); n != NULL; n=n->prev())
 {
  s = ((Triangle *)n->data);
  unlinkTriangle(s);
 }

 removeUnlinkedElements();
}


//////// Next region's boundary vertex /////////////

Vertex *Basic_TMesh::nextVertexOnRegionBoundary(Vertex *sv) const
{
 Triangle *lt, *rt;
 Edge *e;
 List *ve = sv->VE();
 Node *n;

 FOREACHVEEDGE(ve, e, n)
 {
  lt = e->leftTriangle(sv);
  rt = e->rightTriangle(sv);
  if (lt != NULL && IS_VISITED(lt) && (rt == NULL || !IS_VISITED(rt)))
   {delete(ve); return e->oppositeVertex(sv);}
 }
 delete(ve);

 return NULL;
}


//// This method returns a list containing an edge of the region's boundary ////
//// as its first element, and all the internal vertices as the remaining   ////

List *Basic_TMesh::getRegionInternalVertices(List *reg)
{
 List *iVertices = new List;
 List *outList = new List;
 Edge *bEdge = NULL;
 Triangle *s, *t;
 Node *n;
 Vertex *v1, *v2, *v3;

 FOREACHVTTRIANGLE(reg, t, n) {MARK_VISIT(t); MARK_BIT(t, 3);}

 FOREACHVTTRIANGLE(reg, t, n)
 {
  if (IS_BIT(t,3))
  {
   UNMARK_BIT(t,3);
   if ((s = t->t1()) != NULL && !IS_VISITED(s)) {bEdge = t->e1; MARK_BIT(t->e1->v1, 3); MARK_BIT(t->e1->v2, 3);}
   if ((s = t->t2()) != NULL && !IS_VISITED(s)) {bEdge = t->e2; MARK_BIT(t->e2->v1, 3); MARK_BIT(t->e2->v2, 3);}
   if ((s = t->t3()) != NULL && !IS_VISITED(s)) {bEdge = t->e3; MARK_BIT(t->e3->v1, 3); MARK_BIT(t->e3->v2, 3);}
  }
 }

 FOREACHVTTRIANGLE(reg, s, n)
 {
  v1 = s->v1(); v2 = s->v2(); v3 = s->v3();
  if (!IS_BIT(v1, 3)) {iVertices->appendHead(v1); MARK_BIT(v1, 3);}
  if (!IS_BIT(v2, 3)) {iVertices->appendHead(v2); MARK_BIT(v2, 3);}
  if (!IS_BIT(v3, 3)) {iVertices->appendHead(v3); MARK_BIT(v3, 3);}
 }
 FOREACHVTTRIANGLE(reg, s, n)
 {
  v1 = s->v1(); v2 = s->v2(); v3 = s->v3();
  UNMARK_BIT(v1, 3); UNMARK_BIT(v2, 3); UNMARK_BIT(v3, 3); 
 }

 outList->appendHead(iVertices);
 outList->appendHead(bEdge);

 return outList;
}


// Transforms only the shell indicated by 't0'

void Basic_TMesh::transformShell(Triangle *t0, const Matrix4x4& m)
{
 List todo(t0), st, sv;
 Triangle *t, *nt;
 Vertex *v;
 coord x, y, z, w;

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  st.appendHead(t);
  nt=t->t1(); if (nt != NULL && !IS_VISITED(nt)) {MARK_VISIT(nt); todo.appendHead(nt);}
  nt=t->t2(); if (nt != NULL && !IS_VISITED(nt)) {MARK_VISIT(nt); todo.appendHead(nt);}
  nt=t->t3(); if (nt != NULL && !IS_VISITED(nt)) {MARK_VISIT(nt); todo.appendHead(nt);}
 }

 while (st.numels())
 {
  t = (Triangle *)st.popHead();
  UNMARK_VISIT(t);
  v = t->v1(); if (!IS_VISITED(v)) {MARK_VISIT(v); sv.appendHead(v);}
  v = t->v2(); if (!IS_VISITED(v)) {MARK_VISIT(v); sv.appendHead(v);}
  v = t->v3(); if (!IS_VISITED(v)) {MARK_VISIT(v); sv.appendHead(v);}
 }

 while (sv.numels())
 {
  v = (Vertex *)sv.popHead();
  UNMARK_VISIT(v);
  x = ((*v)*Point(m.matrix[0][0],m.matrix[1][0],m.matrix[2][0]))+m.matrix[3][0];
  y = ((*v)*Point(m.matrix[0][1],m.matrix[1][1],m.matrix[2][1]))+m.matrix[3][1];
  z = ((*v)*Point(m.matrix[0][2],m.matrix[1][2],m.matrix[2][2]))+m.matrix[3][2];
  w = ((*v)*Point(m.matrix[0][3],m.matrix[1][3],m.matrix[2][3]))+m.matrix[3][3];
  v->x = x/w; v->y = y/w; v->z = z/w; 
 }
}


//// Removes all the triangles of the shell containing 't0' ////

void Basic_TMesh::removeShell(Triangle *t0)
{
 List todo(t0);
 Triangle *t, *t1, *t2, *t3;

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

  if (t1 != NULL && !IS_VISITED2(t1)) {MARK_VISIT2(t1); todo.appendHead(t1);}
  if (t2 != NULL && !IS_VISITED2(t2)) {MARK_VISIT2(t2); todo.appendHead(t2);}
  if (t3 != NULL && !IS_VISITED2(t3)) {MARK_VISIT2(t3); todo.appendHead(t3);}

  unlinkTriangle(t);
 }

 removeUnlinkedElements();
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    G L O B A L   O P E R A T I O N S                         //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Tags as sharp all the edges exceeding the given curvature ///////////

void Basic_TMesh::sharpEdgeTagging(const double ta)
{
 Node *n;
 Edge *e;
 FOREACHEDGE(e, n)
  if (e->curvature() > ta) TAG_SHARPEDGE(e);
  else UNTAG_SHARPEDGE(e);
}


//// Unmarks all the elements of the triangulation ////

void Basic_TMesh::unmarkEverything()
{
 Vertex *v;
 Edge *e;
 Triangle *t;
 Node *n;
 FOREACHVERTEX(v, n) v->mask = 0;
 FOREACHEDGE(e, n) e->mask = 0;
 FOREACHTRIANGLE(t, n) t->mask = 0;
}


///// Compute the bounding box and return its max edge /////

coord Basic_TMesh::getBoundingBox(Point& mp, Point& Mp) const
{
 Vertex *v; Node *n;
 Mp.x = -DBL_MAX, mp.x = DBL_MAX;
 Mp.y = -DBL_MAX, mp.y = DBL_MAX;
 Mp.z = -DBL_MAX, mp.z = DBL_MAX;
 FOREACHVERTEX(v, n)
 {
  if (v->x < mp.x) mp.x = v->x;
  if (v->x > Mp.x) Mp.x = v->x;
  if (v->y < mp.y) mp.y = v->y;
  if (v->y > Mp.y) Mp.y = v->y;
  if (v->z < mp.z) mp.z = v->z;
  if (v->z > Mp.z) Mp.z = v->z;
 }

 return MAX(Mp.x-mp.x,MAX(Mp.y-mp.y,Mp.z-mp.z));
}


///// Compute the approximate bounding ball radius /////

double Basic_TMesh::getBoundingBallRadius() const
{
 Vertex *v; Node *n;
 Point tc, mp, Mp;
 double tb, bsr = TMESH_TO_DOUBLE(getBoundingBox(mp, Mp)) / 2;
 Point bsc = (Mp+mp)/2;

 FOREACHVERTEX(v, n)
  if ((tb = ((*v)-bsc).length()) > bsr)
  {
   tc = ((*v)-bsc); tc.normalize();
   tb = ((tb-bsr)/2);
   bsc = bsc+(tc*tb);
   bsr += tb;
  }

 return bsr;
}


////// Returns the surface area of the mesh //////

double Basic_TMesh::area() const
{
 Triangle *t;
 Node *n;
 double a=0.0;
 FOREACHTRIANGLE(t, n) a += t->area();

 return a;
}


////// Returns the volume of the mesh //////

double Basic_TMesh::volume() const
{
 Triangle *t;
 Node *n;
 double v=0.0;
 FOREACHTRIANGLE(t, n)
	 v += TMESH_TO_DOUBLE((t->getCenter()*t->getNormal()))*t->area();

 return v/3;
}


///// Places the mesh into the unit cube by translating and resizing /////
///// so that all the coordinates are between 0 and mc (default =1). /////

void Basic_TMesh::normalize(coord mc)
{
	Vertex *v;
	Node *n;
	Point mp, Mp;
	coord mel = getBoundingBox(mp, Mp) / mc;
	FOREACHVERTEX(v, n) v->setValue(((*v) - mp) / mel);	// Shift and normalize
}

void Basic_TMesh::quantize(int nc)
{
	Vertex *v;
	Node *n;
	normalize(nc);
	FOREACHVERTEX(v, n)
	{
		v->x = coord(TMESH_TO_INT(v->x));
		v->y = coord(TMESH_TO_INT(v->y));
		v->z = coord(TMESH_TO_INT(v->z));
	}
}


/////// Transforms all the vertices using the 4x4 matrix 'm' ///////////

void Basic_TMesh::transform(const Matrix4x4& m)
{
 Node *n;
 Vertex *v;
 coord x,y,z,w;

 FOREACHVERTEX(v, n)
 {
  x = ((*v)*Point(m.matrix[0][0],m.matrix[1][0],m.matrix[2][0]))+m.matrix[3][0];
  y = ((*v)*Point(m.matrix[0][1],m.matrix[1][1],m.matrix[2][1]))+m.matrix[3][1];
  z = ((*v)*Point(m.matrix[0][2],m.matrix[1][2],m.matrix[2][2]))+m.matrix[3][2];
  w = ((*v)*Point(m.matrix[0][3],m.matrix[1][3],m.matrix[2][3]))+m.matrix[3][3];
  v->x = x/w; v->y = y/w; v->z = z/w; 
 }
}

// Translate the whole mesh by t_vec
void Basic_TMesh::translate(const Point& t_vec)
{
 Vertex *v;
 Node *i;
 FOREACHVERTEX(v, i) v->setValue((*v)+t_vec);
}

// Return the center of mass of the mesh
Point Basic_TMesh::getCenter() const
{
 Point c;
 Triangle *t;
 Node *i;
 coord av, tvol=0.0;
 FOREACHTRIANGLE(t, i)
 {
  av = t->area();
  tvol += av;
  c += (t->getCenter()*av);
 }

 return c/tvol;
}

// Add noise in the normal direction. Normal displacement is
// bounded by ns% of the bounding ball radius.

void Basic_TMesh::addNormalNoise(double ns)
{
 Vertex *v;
 Node *n;
 Point np;
 int i;
 double noise;
 coord *xyz = (coord *)malloc(sizeof(coord)*V.numels()*3);
 ns *= (getBoundingBallRadius()/100.0);

 i=0; FOREACHVERTEX(v, n)
 {
  noise = ns*(((((double)rand()))-(((double)RAND_MAX)/2.0))/((double)RAND_MAX));
  np = (*v)+((v->getNormal())*noise);
  xyz[i++]=np.x; xyz[i++]=np.y; xyz[i++]=np.z;
 }
 i=0; FOREACHVERTEX(v, n) {v->x=xyz[i++]; v->y=xyz[i++]; v->z=xyz[i++];}

 free(xyz);
}


// Iteratively swap edges to maximize the minimum angle.
// Checks and avoids normal inversion.

bool Basic_TMesh::iterativeEdgeSwaps()
{
 Node *n;
 Edge *e, *f;
 double l;
 int swaps=1, totits=1;
 Point n1, n2, nor;
 List toswap;

 bool selection=0;
 Triangle *t;
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {selection=1; break;}

 FOREACHEDGE(e, n) if (!IS_SHARPEDGE(e) && !e->isOnBoundary())
 {
  MARK_VISIT(e); if ((!selection || (IS_VISITED(e->t1) && IS_VISITED(e->t2)))) toswap.appendTail(e);
 }

 TMesh::begin_progress();
 while (swaps && totits++ < 10)
 {
  swaps = 0; for (n=toswap.head(); n!=NULL; )
  {
   e = (Edge *)n->data;
   if (n==toswap.tail()) {toswap.removeCell(toswap.tail()); n=NULL;}
   else {n=n->next(); toswap.removeCell(n->prev());}
   UNMARK_VISIT(e);
 
    n1 = e->t1->getNormal();
    n2 = e->t2->getNormal();
    nor = n1+n2;
    l = e->delaunayMinAngle();
    if (e->swap())
    {
     if (e->delaunayMinAngle() <= l*1.000001 || nor*e->t1->getNormal() <= 0 || nor*e->t2->getNormal() <= 0) e->swap(1);
     else
     {
      swaps++;
      f = e->t1->nextEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
      f = e->t1->prevEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
      f = e->t2->nextEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
      f = e->t2->prevEdge(e); if (!IS_VISITED(f) && !IS_SHARPEDGE(f) && !f->isOnBoundary()) {MARK_VISIT(f); toswap.appendHead(f);}
     }
    }

  }
  TMesh::report_progress("Swaps: %d      ", swaps);
 }
 TMesh::end_progress();

 FOREACHEDGE(e, n) UNMARK_VISIT(e);

 if (totits >= 10)
 {
  TMesh::warning("Optimization did not converge after 10 iterations! Stopping.\n");
  TMesh::warning("You may try to run the method again.\n");
  return 0;
 }

 return 1;
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    T O P O L O G Y   M A N I P U L A T I O N                 //
//                                                              //
//////////////////////////////////////////////////////////////////


////// Invertes all the triangle normals and edge orientations ///////

void Basic_TMesh::flipNormals()
{
 Node *n;
 Edge *e;
 Triangle *t;

 FOREACHTRIANGLE(t, n) t->invert();
 FOREACHEDGE(e, n) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));		//!< AMF_CHANGE 1.1-2 >
}


//// Marks all the triangles of the shell containing 't' ////

void Basic_TMesh::flipNormals(Triangle *t0)
{
 List todo;
 Triangle *t, *t1, *t2, *t3;

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (!IS_BIT(t,6))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && !IS_BIT(t1,6)) todo.appendHead(t1);
   if (t2 != NULL && !IS_BIT(t2,6)) todo.appendHead(t2);
   if (t3 != NULL && !IS_BIT(t3,6)) todo.appendHead(t3);

   t->invert();
   if (!IS_BIT(t->e1,6)) p_swap((void **)(&(t->e1->v1)), (void **)(&(t->e1->v2)));
   if (!IS_BIT(t->e2,6)) p_swap((void **)(&(t->e2->v1)), (void **)(&(t->e2->v2)));
   if (!IS_BIT(t->e3,6)) p_swap((void **)(&(t->e3->v1)), (void **)(&(t->e3->v2)));
   MARK_BIT(t->e1,6); MARK_BIT(t->e2,6); MARK_BIT(t->e3,6);
   MARK_BIT(t,6);
  }
 }

 todo.appendHead(t0);
 while (todo.numels())
 {
  t = (Triangle *)todo.popHead();
  if (IS_BIT(t,6))
  {
   t1 = t->t1(); t2 = t->t2(); t3 = t->t3();

   if (t1 != NULL && IS_BIT(t1,6)) todo.appendHead(t1);
   if (t2 != NULL && IS_BIT(t2,6)) todo.appendHead(t2);
   if (t3 != NULL && IS_BIT(t3,6)) todo.appendHead(t3);

   UNMARK_BIT(t->e1,6); UNMARK_BIT(t->e2,6); UNMARK_BIT(t->e3,6);
   UNMARK_BIT(t,6);
  }
 }
}


//// Returns the top triangle of the mesh (max. z) ////

Triangle *Basic_TMesh::topTriangle(Triangle *t0)
{
 Node *n;
 Vertex *v, *hv = NULL, *v1, *v2, *v3;
 Edge *e, *fe = NULL;
 coord az, Mz = -DBL_MAX;
 Triangle *t, *t1, *t2, *t3;
 List *ve, todo, tlist, elist, vlist;

 todo.appendHead(t0); MARK_BIT(t0,2);

 while (todo.numels())
 {
  t = (Triangle *)todo.popHead(); tlist.appendHead(t);
  t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
  v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
  if (!IS_VISITED(v1)) {MARK_VISIT(v1); vlist.appendHead(v1);}
  if (!IS_VISITED(v2)) {MARK_VISIT(v2); vlist.appendHead(v2);}
  if (!IS_VISITED(v3)) {MARK_VISIT(v3); vlist.appendHead(v3);}

  if (!IS_VISITED(t->e1)) {MARK_VISIT(t->e1); elist.appendHead(t->e1);}
  if (!IS_VISITED(t->e2)) {MARK_VISIT(t->e2); elist.appendHead(t->e2);}
  if (!IS_VISITED(t->e3)) {MARK_VISIT(t->e3); elist.appendHead(t->e3);}

  if (t1 != NULL && !IS_BIT(t1,2)) {MARK_BIT(t1,2); todo.appendHead(t1);}
  if (t2 != NULL && !IS_BIT(t2,2)) {MARK_BIT(t2,2); todo.appendHead(t2);}
  if (t3 != NULL && !IS_BIT(t3,2)) {MARK_BIT(t3,2); todo.appendHead(t3);}
 }

 ve = new List;

 FOREACHVVVERTEX((&(vlist)), v, n) {UNMARK_VISIT(v); if ((az = v->z) > Mz) {Mz=az; hv = v;}}
 Mz = DBL_MAX;
 FOREACHVEEDGE((&(elist)), e, n) {UNMARK_VISIT(e); if (e->hasVertex(hv) && e->length() != 0) ve->appendHead(e);}
 FOREACHVTTRIANGLE((&(tlist)), t, n) UNMARK_BIT(t, 2);

 FOREACHVEEDGE(ve, e, n)
  if ((az = (hv->z - e->oppositeVertex(hv)->z)/e->length()) < Mz) {Mz=az; fe = e;}
 delete(ve);

 if (fe == NULL) fe = hv->e0;
 if (fe->t1 == NULL || fe->t2 == NULL) return NULL;

 return (FABS(fe->t1->getNormal().z) > FABS(fe->t2->getNormal().z))?(fe->t1):(fe->t2);
}


/////// Computes boundaries and handles ///////////

void Basic_TMesh::eulerUpdate()
{
	Vertex *v, *w;
	Edge *e;
	Triangle *t, *s;
	List triList;
	Node *n;
	n_boundaries = n_shells = n_handles = 0;

	FOREACHTRIANGLE(t, n)	UNMARK_BIT(t, 5);	
	FOREACHVERTEX(v, n)	UNMARK_BIT(v, 5);	

	FOREACHTRIANGLE(t, n) if (!IS_BIT(t, 5))	
	{
		n_shells++;
		triList.appendHead(t);
		MARK_BIT(t, 5);

		while (triList.numels())
		{
			t = (Triangle *)triList.popHead();
			if ((s = t->t1()) != NULL && !IS_BIT(s, 5)) { triList.appendHead(s); MARK_BIT(s, 5); }	
			if ((s = t->t2()) != NULL && !IS_BIT(s, 5)) { triList.appendHead(s); MARK_BIT(s, 5); }	
			if ((s = t->t3()) != NULL && !IS_BIT(s, 5)) { triList.appendHead(s); MARK_BIT(s, 5); }	
		}
	}
	FOREACHTRIANGLE(t, n) UNMARK_BIT(t, 5);	

	bool hasBoundary = false;	

	FOREACHEDGE(e, n) if (e->isOnBoundary()){
		hasBoundary = true;	
		MARK_BIT(e->v1, 5);	
		MARK_BIT(e->v2, 5);	
	}

	if (hasBoundary){			
		FOREACHVERTEX(v, n) if (IS_BIT(v, 5))	
		{
			n_boundaries++;
			for (w = v; IS_BIT(w, 5); w = w->nextOnBoundary()) UNMARK_BIT(w, 5);	
		}
	}

	n_handles = (E.numels() - V.numels() - T.numels() + 2 * n_shells - n_boundaries) / 2;
	d_boundaries = d_handles = d_shells = 0;
}


// Makes the mesh equivalent to a topological disk.
// Edges and vertices are duplicated when necessary.

void Basic_TMesh::openToDisk()
{
 Triangle *t = (Triangle *)T.head()->data;
 Triangle *s;
 List triList, *ve;
 Vertex *v, *w;
 Edge *e, *ne;
 Node *n;
 triList.appendHead(t);
 MARK_BIT(t,3);

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  if ((s = t->t1()) != NULL && !IS_BIT(s,3))
   {triList.appendTail(s); MARK_BIT(s,3); MARK_BIT(t->e1,3);}
  if ((s = t->t2()) != NULL && !IS_BIT(s,3))
   {triList.appendTail(s); MARK_BIT(s,3); MARK_BIT(t->e2,3);}
  if ((s = t->t3()) != NULL && !IS_BIT(s,3))
   {triList.appendTail(s); MARK_BIT(s,3); MARK_BIT(t->e3,3);}
 }
 FOREACHTRIANGLE (t, n) UNMARK_BIT(t, 3);

 FOREACHVERTEX(v, n) v->info = new List;

 FOREACHEDGE(e, n) if (!IS_BIT(e, 3))
 {
  ((List *)e->v1->info)->appendHead(e);
  ((List *)e->v2->info)->appendHead(e);
 }

 FOREACHVERTEX(v, n) if (((List *)v->info)->numels()==1) triList.appendHead(v);
 if (!triList.numels()) TMesh::error("Basic_TMesh::openToDisk: Couldn't find a root.\n");

 while(triList.numels())
 {
  v = (Vertex *)triList.popHead();
  ve = ((List *)v->info);
  if (ve->numels())
  {
   e = (Edge *)(ve->head()->data);
   MARK_BIT(e, 3);
   ve->popHead();
   w = e->oppositeVertex(v);
   ve = ((List *)w->info);
   ve->removeNode(e);
   if (ve->numels() == 1) triList.appendHead(w);
  }
  else
  {
   ve = v->VE();
   e = (Edge *)ve->head()->data; UNMARK_BIT(e, 3); ((List *)v->info)->appendHead(e);
   e = (Edge *)ve->head()->next()->data; UNMARK_BIT(e, 3); ((List *)v->info)->appendHead(e);
   delete(ve);
  }
 }

 FOREACHEDGE(e, n) if (!IS_BIT(e, 3) && !e->isOnBoundary())
 {
  ne = newEdge(e->v1, e->v2);
  ne->t1 = e->t1; e->t1 = NULL; E.appendHead(ne);
  ne->t1->replaceEdge(e, ne);
 }

 FOREACHEDGE(e, n) UNMARK_BIT(e, 3);

 FOREACHVERTEX(v, n) if (v->info) {delete(((List *)v->info)); v->info = NULL;}
 duplicateNonManifoldVertices();
 d_boundaries = d_handles = d_shells = 1;
}


//// Splits the edge 'e' by the point 'p'. ////

Vertex *Basic_TMesh::splitEdge(Edge *e, Point *p, bool copy_mask)
{
 if ((*p)==(*(e->v1))) return e->v1;
 if ((*p)==(*(e->v2))) return e->v2;
 Vertex *v3 = (e->t1 != NULL)?(e->t1->oppositeVertex(e)):(NULL);
 Vertex *v4 = (e->t2 != NULL)?(e->t2->oppositeVertex(e)):(NULL);
 Edge *be1 = (e->t1 != NULL)?(e->t1->nextEdge(e)):(NULL);
 Edge *be4 = (e->t2 != NULL)?(e->t2->prevEdge(e)):(NULL);
 Vertex *v = newVertex(p->x, p->y, p->z);
 Edge *ne = newEdge(v, e->v2);
 Edge *ne1 = (e->t1 != NULL)?(newEdge(v, v3)):(NULL);
 Edge *ne2 = (e->t2 != NULL)?(newEdge(v, v4)):(NULL);
 Triangle *nt1 = (e->t1 != NULL)?(newTriangle(ne1, ne,be1)):(NULL);
 Triangle *nt2 = (e->t2 != NULL)?(newTriangle(ne, ne2,be4)):(NULL);

 ne->t1 = nt1; ne->t2 = nt2;
 if (ne1 != NULL) {ne1->t1 = e->t1; ne1->t2 = nt1;}
 if (ne2 != NULL) {ne2->t1 = nt2; ne2->t2 = e->t2;}
 if (be1 != NULL) be1->replaceTriangle(e->t1, nt1);
 if (be4 != NULL) be4->replaceTriangle(e->t2, nt2);
 e->v2->e0 = (be1 != NULL)?(be1):(be4);
 e->v2 = v;
 v->e0 = e;
 if (e->t1 != NULL) e->t1->replaceEdge(be1, ne1);
 if (e->t2 != NULL) e->t2->replaceEdge(be4, ne2);

 if (copy_mask)
 {
  ne->mask = e->mask;
  if (nt1 != NULL) nt1->mask = e->t1->mask;
  if (nt2 != NULL) nt2->mask = e->t2->mask;
 }

 V.appendHead(v);
 E.appendHead(ne);
 if (ne1 != NULL) E.appendHead(ne1);
 if (ne2 != NULL) E.appendHead(ne2);
 if (nt1 != NULL) T.appendHead(nt1);
 if (nt2 != NULL) T.appendHead(nt2);

 return v;
}


//// Splits the trianlge 't' by the point 'p'. If 'corr' is TRUE ////
//// the method does not perform the split if this causes the   ////
//// creation of a triangle with a vertex angle < MIN_ANG_TRI.  ////

Vertex *Basic_TMesh::splitTriangle(Triangle *t, Point *p, bool copy_mask)
{
 Vertex *v1 = t->v1();
 Vertex *v2 = t->v2();
 Vertex *v3 = t->v3();

 Vertex *v = newVertex(p->x, p->y, p->z);
 Edge *ne1 = newEdge(v, v1);
 Edge *ne2 = newEdge(v, v2);
 Edge *ne3 = newEdge(v, v3);
 Triangle *nt1 = newTriangle(ne2, t->e3,ne3);
 Triangle *nt2 = newTriangle(ne3, t->e1,ne1);
 t->e3->replaceTriangle(t, nt1);
 t->e1->replaceTriangle(t, nt2);
 t->replaceEdge(t->e3, ne2);
 t->replaceEdge(t->e1, ne1);
 ne1->t1 = t; ne1->t2 = nt2;
 ne2->t1 = nt1; ne2->t2 = t;
 ne3->t1 = nt2; ne3->t2 = nt1;
 v->e0 = ne1;

 V.appendHead(v);
 E.appendHead(ne1);
 E.appendHead(ne2);
 E.appendHead(ne3);
 T.appendHead(nt1);
 T.appendHead(nt2);

 if (copy_mask)
 {
	 nt1->mask = t->mask;
	 nt2->mask = t->mask;
 }

 return v;
}


//////////////////////////////////////////////////////////////////
//                                                              //
//    D E B U G   AND WORK-IN-PROGRESS                          //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Prints general information ///////////

void Basic_TMesh::printReport()
{
 eulerUpdate();

 TMesh::info("*** Basic_TMesh Report ***\n");
 TMesh::info("V: %d\n",V.numels());
 TMesh::info("E: %d\n",E.numels());
 TMesh::info("T: %d\n",T.numels());

 TMesh::info("Boundary: %d components.\n",boundaries());
 TMesh::info("Handles: %d.\n",handles());
 TMesh::info("Shells: %d.\n",shells());
}


// Unless there are bugs, the following should be exact ...

bool Basic_TMesh::isInnerPoint(Point& p) const
{
	if (T.numels() == 0) return false; // An ampty mesh does not enclose anything ...

	// We assume that the mesh is correctly oriented.
	Node *n;
	Vertex *v1, *v2, *v3;
	Triangle *t;

	// Ray casting along positive X direction
	Point p2(1, 0, 0);
	p2 += p;
	coord ad, min_distance = DBL_MAX;
	Point ip;
	Triangle *closest_triangle = NULL;
	Edge *e, *closest_edge = NULL;
	Vertex *closest_vertex = NULL;

	double dcp[2]; dcp[0] = TMESH_TO_DOUBLE(p.y); dcp[1] = TMESH_TO_DOUBLE(p.z);
	double dcv[6];
	coord o1, o2, o3;
	FOREACHTRIANGLE(t, n)
	{
		v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
		if (((v1->y > p.y && v2->y > p.y) || (v1->y < p.y && v2->y < p.y)) && ((v1->y > p.y && v3->y > p.y) || (v1->y < p.y && v3->y < p.y)))
			continue;
		if (((v1->z > p.z && v2->z > p.z) || (v1->z < p.z && v2->z < p.z)) && ((v1->z > p.z && v3->z > p.z) || (v1->z < p.z && v3->z < p.z)))
			continue;
#ifdef USE_HYBRID_KERNEL
		if (coord::use_rationals)
		{
			o1 = coord::orient2D(p.y, p.z, v1->y, v1->z, v2->y, v2->z);
			o2 = coord::orient2D(p.y, p.z, v2->y, v2->z, v3->y, v3->z);
			o3 = coord::orient2D(p.y, p.z, v3->y, v3->z, v1->y, v1->z);
		}
		else
		{
#endif
			dcv[0] = TMESH_TO_DOUBLE(v1->y); dcv[1] = TMESH_TO_DOUBLE(v1->z); dcv[2] = TMESH_TO_DOUBLE(v2->y); dcv[3] = TMESH_TO_DOUBLE(v2->z); dcv[4] = TMESH_TO_DOUBLE(v3->y); dcv[5] = TMESH_TO_DOUBLE(v3->z);
			o1 = orient2d(dcp, dcv, dcv + 2);
			o2 = orient2d(dcp, dcv + 2, dcv + 4);
			o3 = orient2d(dcp, dcv + 4, dcv);
#ifdef USE_HYBRID_KERNEL
		}
#endif

		if (o1 == 0 && o2 == 0 && o3 == 0) continue; // Degenerate triangle. Skip.
		else if (o1 == 0 && o2 == 0)
		{
			if ((ad = v2->x - p.x) == 0) return false; // Point is on surface
			else if (ad > 0 && ad<min_distance) { closest_vertex = v2; closest_edge = NULL; closest_triangle = NULL; min_distance = ad; }
		}
		else if (o1 == 0 && o3 == 0)
		{
			if ((ad = v1->x - p.x) == 0) return false; // Point is on surface
			else if (ad > 0 && ad<min_distance) { closest_vertex = v1; closest_edge = NULL; closest_triangle = NULL; min_distance = ad; }
		}
		else if (o2 == 0 && o3 == 0)
		{
			if ((ad = v3->x - p.x) == 0) return false; // Point is on surface
			else if (ad > 0 && ad<min_distance) { closest_vertex = v3; closest_edge = NULL; closest_triangle = NULL; min_distance = ad; }
		}
		else if (o1 == 0 || o2 == 0 || o3 == 0)
		{
			e = (o1 == 0) ? (t->e2) : ((o2 == 0) ? (t->e3) : (t->e1));
			if ((p.y < e->v1->y && p.y < e->v2->y) || (p.y > e->v1->y && p.y > e->v2->y) ||
				(p.z < e->v1->z && p.z < e->v2->z) || (p.z > e->v1->z && p.z > e->v2->z)) continue;
			ip = Point::lineLineIntersection(p, p2, *e->v1, *e->v2);
			ad = ip.x - p.x;
			if (ad == 0) return false; // Point is on surface
			else if (ad > 0 && ad<min_distance) { closest_vertex = NULL; closest_edge = e; closest_triangle = NULL; min_distance = ad; }
		}
		else if ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0))
		{
			ip = Point::linePlaneIntersection(p, p2, *v1, *v2, *v3);
			ad = ip.x - p.x;
			if (ad == 0) return false; // Point is on surface
			else if (ad > 0 && ad<min_distance) { closest_vertex = NULL; closest_edge = NULL; closest_triangle = t; min_distance = ad; }
		}
	}

	if (closest_vertex != NULL)
	{
		List *ve = closest_vertex->VE();
		coord ad, mind = DBL_MAX;
		Point a = p - (*closest_vertex);
		FOREACHVEEDGE(ve, e, n)
		{
			Vertex *ov1 = e->oppositeVertex(closest_vertex);
			Point b = (*ov1) - (*closest_vertex);
			ad = (((a&b)*(a&b)) / ((b*b)*(a*a))) - 1;
			if (a*b < 0) ad = -ad;
			if (ad < mind) { mind = ad; closest_edge = e; }
		}
		delete ve;
	}

	// If cp is in the interior of an edge, select one of its incident triangles
	if (closest_edge != NULL)
	{
		if (closest_edge->isOnBoundary()) return false;
		Vertex *ov1 = closest_edge->t1->oppositeVertex(closest_edge);
		Vertex *ov2 = closest_edge->t2->oppositeVertex(closest_edge);
		coord o1 = p.exactOrientation(closest_edge->v1, closest_edge->v2, ov1);
		coord o2 = ov2->exactOrientation(closest_edge->v1, closest_edge->v2, ov1);
		closest_triangle = ((o1 >= 0 && o2 >= 0) || (o1 <= 0 && o2 <= 0)) ? (closest_edge->t2) : (closest_edge->t1);
	}

	// If closest point is in the interior of a triangle, just check orientation
	if (closest_triangle != NULL)
	{
		return (closest_triangle->getVector().x > 0);
	}

	return false;

	//if (!singular_case)
	//{
	//	int parity = 0;
	//	Point p2(1, 0, 0);
	//	p2 += p;
	//	FOREACHTRIANGLE(t, n)
	//	{
	//		v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
	//		if ((IS_BIT(v1, 5) == IS_BIT(v2, 5) && IS_BIT(v1, 5) == IS_BIT(v3, 5))) continue;
	//		if ((IS_BIT(v1, 6) == IS_BIT(v2, 6) && IS_BIT(v1, 6) == IS_BIT(v3, 6))) continue;
	//		coord o1 = coord::orient2D(p.y, p.z, v1->y, v1->z, v2->y, v2->z);
	//		coord o2 = coord::orient2D(p.y, p.z, v2->y, v2->z, v3->y, v3->z);
	//		coord o3 = coord::orient2D(p.y, p.z, v3->y, v3->z, v1->y, v1->z);
	//		if (o1 == 0 || o2 == 0 || o3 == 0) { singular_case = true; break; }
	//		if ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0))
	//		{
	//			o1 = p.exactOrientation(v1, v2, v3);
	//			if (o1 == 0) return false; // Point is on triangle
	//			Point ip = Point::linePlaneIntersection(p, p2, v1, v2, v3);
	//			if (ip.x>p.x)
	//			{
	//				if (o1 < 0) parity++; else if (o1 > 0) parity--;
	//			}
	//		}
	//	}
	//	if (!singular_case) return (parity == 1);
	//}

	// In case of singularity, revert to the following (slower) method
	// Pick the closest triangle
	//coord ad, min_d = DBL_MAX;
	//Triangle *gt;
	//FOREACHTRIANGLE(t, n) if ((ad = t->pointTriangleSquaredDistance(&p)) < min_d) { gt = t; min_d = ad; }

	//Edge *closest_edge = NULL;
	//Vertex *closest_vertex = NULL;
	//gt->pointTriangleSquaredDistance(&p, &closest_edge, &closest_vertex);

	//// If closest point is in the interior of gt, just check orientation
	//if (closest_edge == NULL) return (p.exactOrientation(gt->v1(), gt->v2(), gt->v3())<0);

	//// If cp is in the interior of an edge, check edge convexity (e->getConvexity())
	//if (closest_vertex == NULL) return (closest_edge->getConvexity() < 0);

	//// If cp is a vertex, find the triangle in VT whose normal is mostly aligned with (cp-p) and check orientation
	//Point normal = closest_vertex->getNormal(); // NOT ROBUST !!! This is an extremely particular case, but ...
	//return (((*closest_vertex-p)*normal)>0);
}

} //namespace T_MESH
