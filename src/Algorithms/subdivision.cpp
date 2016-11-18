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

#include "tmesh.h"

using namespace T_MESH;

//////////////////// Loop's Subdivision Scheme //////////////////

class loopSplit
{
 public:
 Edge *e;
 Point p;

 loopSplit(Edge *_e, int md)
 {
  e = _e;
  if (e->t1 == NULL || e->t2 == NULL || md) p = ((*e->v1)+(*e->v2))/2.0;
  else
  {
   Vertex *ov1 = e->t1->oppositeVertex(e);
   Vertex *ov2 = e->t2->oppositeVertex(e);
   p = ((((*e->v1)+(*e->v2))*3.0)+((*ov1)+(*ov2)))/8.0;
  }
 }
};

double subsurfbeta_loop(int k)
{
 double beta = (cos((2.0*M_PI)/((double)k))/4.0)+(3.0/8.0);
 return ((5.0/8.0)-(beta*beta))/((double)k);
}


void loopRelaxOriginal(Vertex *v)
{
 Node *n;
 Edge *e;
 List *ve;
 Point np;
 int k;
 double beta;

 if (v->isOnBoundary())
 {
  ve = v->VE();
  np = (*(((Edge *)ve->head()->data)->oppositeVertex(v)));
  np += (*(((Edge *)ve->tail()->data)->oppositeVertex(v)));
  np = (((*v)*6.0)+np)/8.0;
  delete(ve);
 }
 else
 {
  ve = v->VE();
  k = ve->numels();
  beta = subsurfbeta_loop(k);
  FOREACHVEEDGE(ve, e, n) np += (*(e->oppositeVertex(v)));
  np = ((*v)*(1.0-k*beta))+(np*beta);
  delete(ve);
 }

 v->info = new Point(&np);
}


//// Performs one subdivision step using Loop's scheme     ////
//// If 'md' is set the method does not alter the geometry ////

void Basic_TMesh::loopSubdivision(bool midpoint)
{
 List nvs;
 Edge *e, *e1, *e2;
 Node *n;
 Vertex *v;
 Triangle *t;
 loopSplit *ls;
 int k, detected_sharp=0, is_selection=0;
 if (!midpoint) FOREACHVERTEX(v, n) loopRelaxOriginal(v);

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {is_selection=1; break;}
 if (is_selection) {FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) {MARK_BIT(t->e1, 3); MARK_BIT(t->e2, 3); MARK_BIT(t->e3, 3);}}
 else FOREACHEDGE(e, n) MARK_BIT(e, 3);

 FOREACHEDGE(e, n) if (IS_BIT(e, 3))
 {
  MARK_BIT(e->v1, 3); MARK_BIT(e->v2, 3);
  if (!midpoint && IS_SHARPEDGE(e)) detected_sharp = 1;
  nvs.appendHead(new loopSplit(e, midpoint));
 }

 if (detected_sharp)
  TMesh::warning("loopSubdivision: Crease-preservation is not supported.\n");

 FOREACHNODE(nvs, n)
 {
  ls = ((loopSplit *)n->data);
  k = ls->e->isOnBoundary();
  v = splitEdge(ls->e, &(ls->p));
  e1 = (Edge *)E.head()->data;
  e2 = (Edge *)E.head()->next()->data;
  MARK_VISIT2(v); MARK_VISIT2(e1); if (!k) MARK_VISIT2(e2);

  if (ls->e->t2)
  {
   if (IS_VISITED(ls->e->t2)) {t=(Triangle *)T.head()->data; MARK_VISIT(t);}
   if (ls->e->t1 && IS_VISITED(ls->e->t1)) {t=(Triangle *)T.head()->next()->data; MARK_VISIT(t);}
  }
  else if (IS_VISITED(ls->e->t1)) {t=(Triangle *)T.head()->data; MARK_VISIT(t);}

  if (IS_SHARPEDGE(ls->e))
  {
   if (k) TAG_SHARPEDGE(e2);
   else TAG_SHARPEDGE((Edge *)E.head()->next()->next()->data);
  }
 }

 nvs.freeNodes();

 FOREACHEDGE(e, n)
  if (IS_VISITED2(e))
  {
   UNMARK_VISIT2(e);
   if ((IS_VISITED2(e->v1) && !IS_VISITED2(e->v2)) ||
       (!IS_VISITED2(e->v1) && IS_VISITED2(e->v2)))
    if (e->swap(1) && (!IS_VISITED2(e->v1) || !IS_VISITED2(e->v2))) e->swap(1);
  }

 FOREACHVERTEX(v, n) if (!IS_VISITED2(v) && !midpoint && IS_BIT(v, 3))
 {
  v->setValue((Point *)v->info);
  delete((Point *)v->info);
  v->info = NULL;
 }
 FOREACHVERTEX(v, n) {UNMARK_VISIT2(v); UNMARK_BIT(v, 3);}

 if (detected_sharp)
  TMesh::warning("loopSubdivision: Tagged sharp edges have been smoothed.\n");

 FOREACHEDGE(e, n) UNMARK_BIT(e, 3);
}
