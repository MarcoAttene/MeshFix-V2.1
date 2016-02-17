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

namespace T_MESH
{

//////////////////////////////////////////////////////////////////
//                                                              //
//    T R I A N G U L A T I O N   M E T H O D S                 //
//                                                              //
//////////////////////////////////////////////////////////////////


/////// Computes the center of mass of the hole and star-patches ////////

int Basic_TMesh::StarTriangulateHole(Edge *e)
{
 if (!e->isOnBoundary()) return 0;

 List bvs;
 Node *n;
 Edge *e1, *e2, *e3;
 Point np;
 Vertex *v, *nv, *v1, *v2;
 int nt=0;

 v = e->v1;

 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
  //MARK_BIT(v,5);		// < AMF_ADD - since IMAT-STL 2.4-1 >
 } while (v != e->v1);

 FOREACHVVVERTEX((&(bvs)), v, n) np = np+(*v);

 np = np/bvs.numels();
 nv = newVertex(&np);
 V.appendHead(nv);

 v1 = ((Vertex *)bvs.head()->data);
 Edge *ep = v1->e0;						// AMF_ADD - since IMAT-STL 2.4-2
 e1 = CreateEdge(nv, v1);
 v1->e0 = ep;							// AMF_ADD - since IMAT-STL 2.4-2
 for (n=bvs.head()->next(); n!=NULL; n=n->next())
 {
  v2 = ((Vertex *)n->data);
  e2 = CreateEdge(nv, v2);
  e3 = v1->getEdge(v2);
  CreateTriangle(e1, e2, e3);
  nt++;
  v1 = v2;
  e1 = e2;
 }
 v2 = ((Vertex *)bvs.head()->data);
 e2 = nv->getEdge(v2);
 e3 = v1->getEdge(v2);
 CreateTriangle(e1, e2, e3);
 nt++;

 return nt;
}


///// Patch holes using 2D Delaunay triangulation on the plane 'nor' /////

int Basic_TMesh::TriangulateHole(Edge *e, Point *nor)
{
 if (!e->isOnBoundary()) return 0;

 List bvs;
 Node *n, *gn = NULL;
 Edge *e1, *e2;
 Vertex *v, *v1, *v2;
 double ang, gang;
 int nt = 0;

 v = e->v1;
 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);

 while (bvs.numels() > 2)
 {
  gang = DBL_MAX;
  FOREACHVVVERTEX((&(bvs)), v, n)
   if (!IS_BIT(v, 5) && v->e0 && (ang = v->getAngleOnAveragePlane(nor)) < gang)
    {gang = ang; gn = n;}
  if (gang == DBL_MAX)
  {
   TMesh::warning("TriangulateHole: Can't complete the triangulation.\n");
   FOREACHVVVERTEX((&(bvs)), v, n) UNMARK_BIT(v, 5);
   return 0;
  }
  v = ((Vertex *)gn->data);
  v1 = (Vertex *)((gn->next() != NULL)?(gn->next()):(bvs.head()))->data;
  v2 = (Vertex *)((gn->prev() != NULL)?(gn->prev()):(bvs.tail()))->data;
  e1 = v->getEdge(v1);
  e2 = v->getEdge(v2);
  if (!EulerEdgeTriangle(e1, e2)) MARK_BIT(v, 5);
  else { bvs.removeCell(gn); UNMARK_BIT(v1, 5); UNMARK_BIT(v2, 5); nt++; }
 }

 int i, skips;
 do
 {
  skips = 0;
  for (n=E.head(), i=2*nt*nt; i<nt; n=n->next(), i--)
  {
   e = ((Edge *)n->data);
   ang = e->delaunayMinAngle();
   if (e->swap())
    {if (e->delaunayMinAngle() <= ang) e->swap(1); else skips++;}
  }
  if (i < 0) {TMesh::warning("Optimization is taking too long. I give up.\n"); break;}
 } while (skips);

 return nt;
}


///// Triangulates a hole using the additional vertices in 'vl' /////

int Basic_TMesh::TriangulateHole(Edge *e, List *vl)
{
 if (!e->isOnBoundary()) return 0;

 List bvs, ovbs, nedg;
 Node *n, *gn = NULL;
 Edge *e1, *e2;
 Vertex *v, *v1, *v2;
 double ang, gang;
 int nt = 0, neb;

 v = e->v1;

 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);
 ovbs.appendList(&bvs);

 while (bvs.numels() > 2) // While there are more than two boundary vertices
 {
  gang = DBL_MAX;
  FOREACHVVVERTEX((&(bvs)), v, n)
   if (!IS_BIT(v, 5) && v->e0 && (ang = v->getAngleForTriangulation()) < gang)
    {gang = ang; gn = n;}
  if (gang == DBL_MAX)
  {
   TMesh::warning("TriangulateHole: Can't complete the triangulation.\n");
   FOREACHVVVERTEX((&(bvs)), v, n) UNMARK_BIT(v, 5);
   return 0;
  }
  v = ((Vertex *)gn->data);
  v1 = (Vertex *)((gn->next() != NULL)?(gn->next()):(bvs.head()))->data;
  v2 = (Vertex *)((gn->prev() != NULL)?(gn->prev()):(bvs.tail()))->data;
  e1 = v->getEdge(v1);
  e2 = v->getEdge(v2);
  neb = E.numels();
  if (!EulerEdgeTriangle(e1, e2)) MARK_BIT(v, 5);
  else
  {
   bvs.removeCell(gn);
   UNMARK_BIT(v1, 5);
   UNMARK_BIT(v2, 5);
   nt++;
   if (E.numels() > neb) nedg.appendHead(E.head()->data);
  }
 }

// if (nt < 2) return nt;

 // Calcolo una normale per il buco come media dei nuovi triangoli
 int i;
 Point nor;

 for (i=0, n=T.head(); i<nt; i++, n=n->next()) nor = nor+((Triangle *)n->data)->getNormal();
 if (nor.isNull())
 {
  TMesh::warning("TriangulateHole: Unable to compute an average normal. Can't optimize.\n");
  return nt;
 }
 nor.normalize();
 
 // Memorizzo da qualche parte la posizione originale dei vertici
 // e Proietto il boundary sul piano con normale quella appena calcolata
 //coord *ovps = (coord *)malloc((ovbs.numels())*3*sizeof(coord));		// AMF_CHANGE - since IMAT-STL 2.4-2
 coord *ovps = new coord[3*(ovbs.numels())];							// AMF_CHANGE - since IMAT-STL 2.4-2
 int j = 0;
 FOREACHVVVERTEX((&(ovbs)), v, n)
 {
  ovps[j++] = v->x; 
  ovps[j++] = v->y; 
  ovps[j++] = v->z;
  v->project(&nor);
 }

 // Proietto i punti interni sul piano d'appoggio
 Point *p;
 //coord *ovpsi = (coord *)malloc((vl->numels())*3*sizeof(coord));		// AMF_CHANGE - since IMAT-STL 2.4-2
 coord *ovpsi = new coord[3*(vl->numels())];							// AMF_CHANGE - since IMAT-STL 2.4-2
 j = 0;
 FOREACHNODE((*vl), n)
 {
  p = ((Point *)n->data);
  ovpsi[j++] = p->x; ovpsi[j++] = p->y; ovpsi[j++] = p->z;
 }
 FOREACHNODE((*vl), n)
 {
  p = ((Point *)n->data);
  p->project(&nor);
 }

 // Ottimizzo secondo Delaunay vincolato al boundary la nuova regione

 int sw;
 do
 {
  sw = 0;
  FOREACHVEEDGE((&(nedg)), e1, n)
  {
   ang = e1->delaunayMinAngle();
   if (e1->swap())
    {if (e1->delaunayMinAngle() <= ang) e1->swap(1); else sw++;}
  }
 } while (sw);

 // Inserisco i punti interni
 int ntt = T.numels()-nt;
 List ivs;

 FOREACHNODE((*vl), n)
 {
  p = ((Point *)n->data);
  ivs.appendTail(watsonInsert(p, &T, T.numels()-ntt));
 }
 nt = (T.numels() - ntt);

 // Riporto i vertici interni al loro posto
 j=0; FOREACHVVVERTEX((&(ivs)), v, n)
 {
  if (v != NULL)
  {
   v->x = ovpsi[j++]; v->y = ovpsi[j++]; v->z = ovpsi[j++];
  } else j+=3;
 }
 delete [] ovpsi; //free(ovpsi);					// AMF_CHANGE - since IMAT-STL 2.4-2

 // Riporto i vertici del boundary al loro posto
 j=0; FOREACHVVVERTEX((&(ovbs)), v, n)
 {
  v->x = ovps[j++]; v->y = ovps[j++]; v->z = ovps[j++];
 }
 delete [] ovps; //free(ovps);						// AMF_CHANGE - since IMAT-STL 2.4-2

 return nt;
}


////// Inserts the point 'p' in the Delaunay triangulation 'tR' with 'nt' triangles ////

Vertex *Basic_TMesh::watsonInsert(Point *p, List *tR, int nt)
{
 Node *n, *m;
 Edge *e;
 Triangle *t;
 List bdr, bdrs, todo, *ve;
 Vertex *v1, *v2, *v3;
 int i;

 for (i=0, n = T.head(); i<nt; n=n->next(), i++)
 {
  t = ((Triangle *)n->data);
  if (t->e1 != NULL && t->inSphere(p))
  {
   v1 = t->v1(); v2 = t->v2(); v3 = t->v3();
   if (!IS_BIT(v1, 5)) bdr.appendHead(v1);
   if (!IS_BIT(v2, 5)) bdr.appendHead(v2);
   if (!IS_BIT(v3, 5)) bdr.appendHead(v3);
   MARK_BIT(v1, 5); MARK_BIT(v2, 5); MARK_BIT(v3, 5);
   MARK_BIT(t, 6);
   todo.appendHead(t);
  }
 }
 if (bdr.numels() == 0) return NULL;

 FOREACHVVVERTEX((&(bdr)), v1, n)
 {
  ve = v1->VE();
  FOREACHVEEDGE(ve, e, m) if (!IS_BIT(e->t1, 6) || !IS_BIT(e->t2, 6)) v1->e0 = e;
  delete(ve);
 }

 while (todo.numels())
 {
  t = ((Triangle *)todo.head()->data);
  todo.removeCell(todo.head());
  unlinkTriangleNoManifold(t);
 }

 Node *tmp;
 for (i=0, n = T.head(); i<nt; i++)
 {
  t = ((Triangle *)n->data);
  if (t->e1 == NULL) {tmp = n; n=n->next(); T.freeCell(tmp);}
  else n=n->next();
 }

 for (n = bdr.head(); n!=NULL;)
 {
  v1 = ((Vertex *)n->data);
  if (v1->e0 == NULL) {tmp = n; n=n->next(); bdr.removeCell(tmp);}
  else n=n->next();
 }

 v1 = v2 = ((Vertex *)bdr.head()->data);
 do
 {
  bdrs.appendHead(v1);
  v1 = v1->nextOnBoundary();
 } while (v1 != v2);

 Vertex *v = newVertex(p->x, p->y, p->z);
 V.appendHead(v);

 v1 = ((Vertex *)bdrs.head()->data);
 v->e0 = e = newEdge(v, v1);
 UNMARK_BIT(v1, 5);
 E.appendHead(e);

 for (n = bdrs.head()->next(); n!=NULL; n=n->next())
 {
  v1 = ((Vertex *)n->data);
  UNMARK_BIT(v1, 5);
  v2 = ((Vertex *)n->prev()->data);
  e = newEdge(v, v1);
  CreateTriangle(e, v1->getEdge(v2), (Edge *)E.head()->data);
  E.appendHead(e);
 }
 EulerEdgeTriangle(v->e0, (Edge *)E.head()->data);

 return v;
}


//// Removes a vertex and retriangulates its VT ////

int Basic_TMesh::retriangulateVT(Vertex *v)
{
 Point nor;
 Edge *e, *e0 = v->e0->t1->oppositeEdge(v);
 List *vt = v->VT();
 List oe;
 Triangle *t;
 Node *m, *n;
 int i, nt;

 FOREACHVTTRIANGLE(vt, t, m)
 {
  e = t->oppositeEdge(v);
  oe.appendTail(t->prevEdge(e));
  oe.appendTail(e);
  oe.appendTail(t->nextEdge(e));
  nor = nor+t->getNormal();
  unlinkTriangle(t);
 }
 removeUnlinkedElements();						// AMF_CHANGE - since IMAT-STL 2.4-2
 nor.normalize();
 nt = TriangulateHole(e0, &nor);

 for (m=T.head(), i=0; i<nt; i++, m=m->next())
 {
  t = ((Triangle *)m->data);
  if (t->overlaps() || t->isExactlyDegenerate()) break;
 }
 if (i<nt)
 {
  TMesh::warning("Re-triangulation failed. Restoring..\n");
  for (m=T.head(), i=0; i<nt; i++, m=m->next())
   unlinkTriangle(((Triangle *)m->data));
  n = oe.head();
  FOREACHVTTRIANGLE(vt, t, m)
  {
   t->e1 = ((Edge *)n->data); n=n->next();
   t->e2 = ((Edge *)n->data); n=n->next();
   t->e3 = ((Edge *)n->data); n=n->next();
   t->e1->v1 = v; t->e1->v2 = (t->e2->t1 == NULL)?(t->e2->v1):(t->e2->v2);
   t->e3->v1 = v; t->e3->v2 = (t->e2->t1 == NULL)?(t->e2->v2):(t->e2->v1);
   ((t->e2->t1 == NULL)?(t->e2->t1):(t->e2->t2)) = t;
   t->e1->t1 = t;
   t->e3->t2 = t;
  }
  v->e0 = ((Triangle *)vt->head()->data)->e1;
 }
 delete(vt);

 return 1;
}


////////// Generic method for patching holes. Heuristic. /////////////
////////// Small angles are patched first, if possible.  /////////////

int Basic_TMesh::TriangulateHole(Edge *e)
{
 if (!e->isOnBoundary()) return 0;

 List bvs;
 Node *n, *gn = NULL;
 Edge *e1, *e2;
 Vertex *v, *v1, *v2;
 double ang, gang;
 int nt = 0;
 Triangle *t;
 v = e->v1;

 t = (e->t1!=NULL)?(e->t1):(e->t2);
 if (t->nextEdge(e)->isOnBoundary() && t->prevEdge(e)->isOnBoundary()) return 0;

 do
 {
  bvs.appendHead(v);
  v = v->nextOnBoundary();
 } while (v != e->v1);

 while (bvs.numels() > 2)
 {
  gang = DBL_MAX;
  FOREACHVVVERTEX((&bvs), v, n)
   if (!IS_BIT(v, 5) && v->e0 && (ang = v->getAngleForTriangulation()) < gang)
    {gang = ang; gn = n;}
  if (gang == DBL_MAX)
  {
   TMesh::warning("TriangulateHole: Can't complete the triangulation.\n");
   FOREACHVVVERTEX((&bvs), v, n) UNMARK_BIT(v, 5);
   int i=0; FOREACHTRIANGLE(t, n) if (i++==nt) break; else unlinkTriangle(t);
   removeUnlinkedElements();
   return 0;
  }
  v = ((Vertex *)gn->data);
  v1 = (Vertex *)((gn->next() != NULL)?(gn->next()):(bvs.head()))->data;
  v2 = (Vertex *)((gn->prev() != NULL)?(gn->prev()):(bvs.tail()))->data;
  e1 = v->getEdge(v1);
  e2 = v->getEdge(v2);
  if ((t=EulerEdgeTriangle(e1,e2))==NULL) MARK_BIT(v, 5);
  else { bvs.removeCell(gn); UNMARK_BIT(v1, 5); UNMARK_BIT(v2, 5); MARK_VISIT(t); nt++; }
 }

 return nt;
}


// Fills the hole identified by 'e' and leaves the new triangle selected.
// 'refine' is for internal vertex insertion.

void Basic_TMesh::FillHole(Edge *e, bool refine)
{
 int i, nt;
 Node *n;
 Triangle *t;
 Vertex *v;

 deselectTriangles();
 FOREACHVERTEX(v, n) UNMARK_BIT(v, 5);
 nt = TriangulateHole(e);
 if (!nt) return;

 i=0; FOREACHTRIANGLE(t, n) if (i++==nt) break; else MARK_VISIT(t);

 if (refine)
  refineSelectedHolePatches((Triangle *)T.head()->data);
}

//// Triangulate Small Boundaries (with less than 'nbe' edges) /////

int Basic_TMesh::fillSmallBoundaries(int nbe, bool refine_patches)
{
 if (nbe == 0) nbe = E.numels();
 Vertex *v,*w;
 Triangle *t;
 Node *n;
 int grd, is_selection=0, tbds = 0, pct = 100;
 List bdrs;

 TMesh::begin_progress();
 TMesh::report_progress("0%% done ");

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {is_selection=1; break;}

 if (is_selection) FOREACHTRIANGLE(t, n) if (!IS_VISITED(t))
  {MARK_BIT(t->v1(), 6); MARK_BIT(t->v2(), 6); MARK_BIT(t->v3(), 6);}

 FOREACHVERTEX(v, n)
 {
  grd = 0;
  if (!IS_BIT(v, 6) && v->isOnBoundary())
  {
   tbds++;
   w = v;
   do
   {
    if (IS_BIT(w, 6)) grd=nbe+1;
	MARK_BIT(w, 6);
    grd++;
    w = w->nextOnBoundary();
   } while (w != v);
   if (grd <= nbe) bdrs.appendHead(w->nextBoundaryEdge());
  }
 }
 FOREACHVERTEX(v, n) { UNMARK_BIT(v, 5); UNMARK_BIT(v, 6); }

 deselectTriangles();

 pct=0; FOREACHNODE(bdrs, n)
 {
  if (TriangulateHole((Edge *)n->data) && refine_patches)
  {
   t = (Triangle *)T.head()->data;
   refineSelectedHolePatches(t);
  }
  TMesh::report_progress("%d%% done ",((++pct)*100)/bdrs.numels());
 }

 grd = bdrs.numels();

 TMesh::end_progress();

 return grd;
}


// Inserts new vertices in the current selection so as
// to reflect the density of the surrounding mesh.
// This method assumes that the selection has no internal vertices.

int Basic_TMesh::refineSelectedHolePatches(Triangle *t0)
{
 Node *n, *m;
 Triangle *t, *t1, *t2;
 Edge *e, *f;
 Vertex *v;
 List *ve, toswap, reg, all_edges, interior_edges, boundary_edges, boundary_vertices, interior_vertices;
 coord sigma, l, sv1, sv2, sv3, dv1, dv2, dv3;
 int swaps, totits, nee, ntb, nnt=-1, pnnt, gits=0;
 const double alpha = sqrt(2.0);
 Point vc;

 if (t0 != NULL)
 {
  if (!IS_VISITED(t0)) TMesh::error("refineSelectedHolePatches: unexpected unselected t0 !");
  UNMARK_VISIT(t0); toswap.appendHead(t0);
  while ((t=(Triangle *)toswap.popHead()) != NULL)
  {
   reg.appendHead(t);
   t1=t->t1(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
   t1=t->t2(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
   t1=t->t3(); if (IS_VISITED(t1)) {UNMARK_VISIT(t1); toswap.appendHead(t1);}
  }
  FOREACHVTTRIANGLE((&reg), t, n) MARK_VISIT(t);
 }
 else FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) reg.appendHead(t);

 printf("%d\n",reg.numels());

 FOREACHVTTRIANGLE((&reg), t, n)
 {
  e = t->e1; if (!IS_BIT(e, 5)) {MARK_BIT(e, 5); all_edges.appendHead(e);} else UNMARK_BIT(e, 5);
  e = t->e2; if (!IS_BIT(e, 5)) {MARK_BIT(e, 5); all_edges.appendHead(e);} else UNMARK_BIT(e, 5);
  e = t->e3; if (!IS_BIT(e, 5)) {MARK_BIT(e, 5); all_edges.appendHead(e);} else UNMARK_BIT(e, 5);
 }

 while (all_edges.numels())
 {
  e = (Edge *)all_edges.popHead();
  if (IS_BIT(e, 5)) { boundary_edges.appendHead(e); UNMARK_BIT(e, 5); }
  else { interior_edges.appendHead(e); MARK_BIT(e, 5); }
 }

 FOREACHVEEDGE((&boundary_edges), e, n)
 {
  v = e->v1; if (!IS_BIT(v, 5)) { MARK_BIT(v, 5); boundary_vertices.appendHead(v);}
  v = e->v2; if (!IS_BIT(v, 5)) { MARK_BIT(v, 5); boundary_vertices.appendHead(v);}
 }

 FOREACHVVVERTEX((&boundary_vertices), v, n) UNMARK_BIT(v, 5);

 // Due to the above definitions, interior edges are BIT

 FOREACHVVVERTEX((&boundary_vertices), v, n)
 {
  ve = v->VE();
  sigma=0; nee=0; FOREACHVEEDGE(ve, e, m) if (!IS_BIT(e, 5)) {nee++; sigma += e->length();}
  sigma /= nee; v->info = new coord(sigma);
  delete(ve);
 }

 FOREACHVEEDGE((&interior_edges), e, n) UNMARK_BIT(e, 5);
 FOREACHVEEDGE((&boundary_edges), e, n) MARK_BIT(e, 6);

 do
 {
  pnnt=nnt;
  nnt=0;
  FOREACHVTTRIANGLE((&reg), t, n)
  {
   vc = t->getCenter();
   sv1 = (*(coord *)t->v1()->info);
   sv2 = (*(coord *)t->v2()->info);
   sv3 = (*(coord *)t->v3()->info);
   sigma = (sv1+sv2+sv3)/3.0;
   dv1 = alpha*(t->v1()->distance(&vc));
   dv2 = alpha*(t->v2()->distance(&vc));
   dv3 = alpha*(t->v3()->distance(&vc));
   if (dv1>sigma && dv1>sv1 && dv2>sigma && dv2>sv2 && dv3>sigma && dv3>sv3)   
   {
    ntb = T.numels();
    v = splitTriangle(t,&vc,1);
    nnt += (T.numels()-ntb);
    if (T.numels() == ntb+2)
    {
     v->info = new coord(sigma);
     interior_vertices.appendHead(v);
     interior_edges.appendHead(v->e0);
     interior_edges.appendHead(v->e0->leftTriangle(v)->prevEdge(v->e0));
     interior_edges.appendHead(v->e0->rightTriangle(v)->nextEdge(v->e0));
     t1 = ((Triangle *)T.head()->data);
     t2 = ((Triangle *)T.head()->next()->data);
     t1->mask = t2->mask = t->mask;
     reg.appendHead(t1); reg.appendHead(t2);
    }
   }
  }

  FOREACHVEEDGE((&interior_edges), e, n) {MARK_BIT(e, 5); toswap.appendHead(e);}
  totits=0; swaps=1;
  while (swaps && totits++ < 10)
  {
   swaps = 0; 
   while ((e=(Edge *)toswap.popHead())!=NULL)
   {
    UNMARK_BIT(e, 5);
    l = e->squaredLength();
    if (e->swap())
    {
     if (e->squaredLength() >= l*0.999999) e->swap(1);
     else
     {
      swaps++;
      toswap.appendTail(e);
	  f = e->t1->nextEdge(e); if (!IS_BIT(f, 5) && !IS_BIT(f, 6)) { MARK_BIT(f, 5); toswap.appendTail(f); }
	  f = e->t1->prevEdge(e); if (!IS_BIT(f, 5) && !IS_BIT(f, 6)) { MARK_BIT(f, 5); toswap.appendTail(f); }
	  f = e->t2->nextEdge(e); if (!IS_BIT(f, 5) && !IS_BIT(f, 6)) { MARK_BIT(f, 5); toswap.appendTail(f); }
	  f = e->t2->prevEdge(e); if (!IS_BIT(f, 5) && !IS_BIT(f, 6)) { MARK_BIT(f, 5); toswap.appendTail(f); }
     }
    }
   }
  }

  if (pnnt==nnt) gits++;
 } while (nnt && gits<10);

 //FOREACHVEEDGE((&boundary_edges), e, n) UNMARK_BIT(e, 6);													
 FOREACHVVVERTEX((&boundary_vertices), v, n) { delete((coord *)v->info); v->info = NULL; MARK_BIT(v, 5);}	
 FOREACHVVVERTEX((&interior_vertices), v, n) { delete((coord *)v->info); v->info = NULL; MARK_BIT(v, 6);}	

 if (gits>=10) {TMesh::warning("Fill holes: Refinement stage failed to converge. Breaking.\n"); return 1;}

 return 0;
}

// Joins the two boundary vertices gv and gw through an edge. A pair of triangles is
// added to properly change the topology of the mesh.
// On success, the return value is the new edge connecting the two vertices.
// NULL is returned on failure.
// Failure occurs if gv and gw are not both on boundary.
// If 'justconnect' is false, the remaining hole is filled with new triangles, unless
// gv and gw are contiguous on the boundary loop (failure).
// If 'justconnect' is true, gv and gw must not belong to the same boundary loop (failure).
// If 'refine', the patching triangles are refined to reproduce neighboring density.

Edge *Basic_TMesh::joinBoundaryLoops(Vertex *gv, Vertex *gw, bool justconnect, bool refine)
{
	Vertex *v, *gvn, *gwn;
	Edge *e, *gve, *gwe;
	Triangle *t;
	Node *n;
	double tl1 = 0.0, tl2 = 0.0, pl1, pl2;

	if (gv == NULL || gw == NULL || !gv->isOnBoundary() || !gw->isOnBoundary()) return NULL;

	FOREACHVERTEX(v, n) UNMARK_VISIT(v);
	deselectTriangles();

	v = gv;
	if (!justconnect)
	{
		do { v = v->nextOnBoundary(); if (v == gw) return NULL; } while (v != gv);
	} else
	{
		gvn = gv->nextOnBoundary(); gwn = gv->prevOnBoundary();
		if (gw == gvn || gw == gwn) return NULL;
		if (gw == gvn->nextOnBoundary())
		{
			t = EulerEdgeTriangle(gvn->prevBoundaryEdge(), gvn->nextBoundaryEdge()); MARK_VISIT(t); return t->oppositeEdge(gvn);
		}
		if (gw == gwn->prevOnBoundary())
		{
			t = EulerEdgeTriangle(gwn->prevBoundaryEdge(), gwn->nextBoundaryEdge()); MARK_VISIT(t); return t->oppositeEdge(gwn);
		}
	}

	gve = gv->prevBoundaryEdge();
	gvn = gve->oppositeVertex(gv);
	gwe = gw->nextBoundaryEdge();
	gwn = gwe->oppositeVertex(gw);

	Edge *je = CreateEdge(gv, gw);
	Edge *je1 = CreateEdge(gv, gwn);
	Edge *je2 = CreateEdge(gwn, gvn);

	t = CreateTriangle(je, gwe, je1); MARK_VISIT(t);
	t = CreateTriangle(je1, je2, gve); MARK_VISIT(t);

	if (justconnect) return je;

	v = gv; do { e = v->nextBoundaryEdge(); v = e->oppositeVertex(v); tl1 += e->length(); } while (v != gv);
	v = gw; do { e = v->nextBoundaryEdge(); v = e->oppositeVertex(v); tl2 += e->length(); } while (v != gw);
	pl1 = tl1; pl2 = tl2;

	double c1, c2;

	e = je;
	while (e->isOnBoundary())
	{
		gv = (e->t2 != NULL) ? (e->v2) : (e->v1); gve = gv->nextBoundaryEdge();
		gw = (e->t1 != NULL) ? (e->v2) : (e->v1); gwe = gw->prevBoundaryEdge();
		c1 = fabs((pl1 - gve->length())*tl2 - pl2*tl1);
		c2 = fabs((pl2 - gwe->length())*tl1 - pl1*tl2);
		if (c1<c2)
		{
			t = EulerEdgeTriangle(e, gve); MARK_VISIT(t);
			pl1 -= gve->length();
			e = t->nextEdge(gve);
		} else
		{
			t = EulerEdgeTriangle(gwe, e); MARK_VISIT(t);
			pl2 -= gwe->length();
			e = t->prevEdge(gwe);
		}
	}

	if (refine) refineSelectedHolePatches();

	return je;
}

} //namespace T_MESH
