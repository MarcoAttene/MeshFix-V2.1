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

#include "detectIntersections.h"
#include <string.h>
#include <stdlib.h>
#include "jqsort.h"

namespace T_MESH
{
	di_cell::di_cell(Basic_TMesh *tin, bool useAll)
	{
		Node *n;
		Vertex *v;
		Triangle *t;
		Mp.x = -DBL_MAX, mp.x = DBL_MAX;
		Mp.y = -DBL_MAX, mp.y = DBL_MAX;
		Mp.z = -DBL_MAX, mp.z = DBL_MAX;
		FOREACHVVVERTEX((&(tin->V)), v, n) if (useAll || IS_BIT(v,5))
		{
			if (v->x < mp.x) mp.x = v->x;
			if (v->x > Mp.x) Mp.x = v->x;
			if (v->y < mp.y) mp.y = v->y;
			if (v->y > Mp.y) Mp.y = v->y;
			if (v->z < mp.z) mp.z = v->z;
			if (v->z > Mp.z) Mp.z = v->z;
		}

		mp -= DI_EPSILON_POINT;
		Mp += DI_EPSILON_POINT;

		FOREACHVTTRIANGLE((&(tin->T)), t, n) if (useAll || IS_VISITED(t)) triangles.appendTail(t);
	}

	bool di_cell::is_triangleBB_in_cell(Triangle *t) const
	{
		Vertex *v1 = t->v1(), *v2 = t->v2(), *v3 = t->v3();
		coord mx = MIN(v1->x, MIN(v2->x, v3->x));
		coord Mx = MAX(v1->x, MAX(v2->x, v3->x));
		coord my = MIN(v1->y, MIN(v2->y, v3->y));
		coord My = MAX(v1->y, MAX(v2->y, v3->y));
		coord mz = MIN(v1->z, MIN(v2->z, v3->z));
		coord Mz = MAX(v1->z, MAX(v2->z, v3->z));

		// Triangle BB is not entirely out of cell
		if (!(Mx<mp.x || mx>Mp.x || My<mp.y || my>Mp.y || Mz<mp.z || mz>Mp.z)) return true;
		else return false;
	}

	di_cell *di_cell::fork()
	{
		Node *n;
		Triangle *t;
		Point e = Mp - mp;
		di_cell *nc = new di_cell;
		char which_coord = 2;

		if (e.x >= e.y && e.x >= e.z) which_coord = 0;
		else if (e.y >= e.x && e.y >= e.z) which_coord = 1;
		nc->mp = mp; nc->Mp = Mp;
		nc->Mp[which_coord] -= (e[which_coord] / 2); mp[which_coord] = nc->Mp[which_coord];

		n = triangles.head();
		while (n != NULL)
		{
			t = (Triangle *)n->data;
			n = n->next();
			if (!is_triangleBB_in_cell(t))
			{
				triangles.moveNodeTo((n != NULL) ? (n->prev()) : triangles.tail(), &(nc->triangles));
			}
			else if (nc->is_triangleBB_in_cell(t)) nc->triangles.appendHead(t);
		}

		return nc;
	}


	// Brute force all-with-all intersection test of the triangles in 'triangles'.
	void di_cell::selectIntersections(bool justproper)
	{
		Triangle *t, *y;
		Node *n, *m;
		List *ts;

		for (n = triangles.head(); n != NULL; n = n->next())
		for (m = n->next(); m != NULL; m = m->next())
		{
			t = (Triangle *)n->data;
			y = (Triangle *)m->data; // For any pair (t,y) of triangles in the cell
			// The same triangle pair can be in different cells. The following avoids redoing the check.
			if (t->info == NULL || y->info == NULL || (((List *)t->info)->containsNode(y) == NULL))
			{
				if (t->intersects(y, justproper))
				{
					MARK_VISIT(t); MARK_VISIT(y);
					ts = ((t->info != NULL) ? ((List *)t->info) : (new List)); ts->appendTail(y); t->info = ts;
					ts = ((y->info != NULL) ? ((List *)y->info) : (new List)); ts->appendTail(t); y->info = ts;
				}
			}
		}
	}


/////////////////////////////////////////////////////////////////////////
//                                                                     ||
////////////////////// Select   Intersections ///////////////////////////
//                                                                     ||
/////////////////////////////////////////////////////////////////////////

int Basic_TMesh::selectIntersectingTriangles(UINT16 tris_per_cell, bool justproper)
{
 Triangle *t;
 Vertex *v;
 Node *n;
 bool isSelection=0;
 List *selT = new List, *selV = new List;
 
 TMesh::begin_progress();
 TMesh::report_progress(NULL);
 
 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
 {
  isSelection=1;
  selT->appendTail(t);
  v=t->v1(); if (!IS_BIT(v,5)) {MARK_BIT(v,5); selV->appendTail(v);}
  v=t->v2(); if (!IS_BIT(v,5)) {MARK_BIT(v,5); selV->appendTail(v);}
  v=t->v3(); if (!IS_BIT(v,5)) {MARK_BIT(v,5); selV->appendTail(v);}
 }
 TMesh::report_progress(NULL);
 
 if (!isSelection) {delete(selT); delete(selV); selT=&T; selV=&V;}

 di_cell *c2, *c = new di_cell(this, !isSelection);
 List cells, todo(c);
 int i=0;

 while ((c = (di_cell *)todo.popHead()) != NULL)
 {
  if (i>DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= tris_per_cell) cells.appendHead(c);
  else
  {
   if (!(i % 1000)) TMesh::report_progress(NULL);
   i++;
   c2 = c->fork();
   todo.appendTail(c);
   todo.appendTail(c2);
  }
 }

 // Deselect everything and select only intersecting triangles
 deselectTriangles();
 FOREACHTRIANGLE(t, n) t->info = NULL;
 i=0; FOREACHNODE(cells, n)
 {
  (((di_cell *)n->data)->selectIntersections(justproper));
  if (!(i % 100)) TMesh::report_progress("%d %% done   ", ((i)* 100) / cells.numels());
  i++;
 }
 TMesh::end_progress();

 // Dispose memory allocated for cells
 FOREACHVTTRIANGLE(selT, t, n) { if (t->info!=NULL) delete((List *)t->info); t->info = NULL; }
 while (cells.numels()) delete((di_cell *)cells.popHead());

 // Count selected triangles for final report and delete stored normals
 int its=0;
 FOREACHVTTRIANGLE(selT, t, n) { if (IS_VISITED(t)) its++;}

 if (its) TMesh::info("%d intersecting triangles have been selected.\n",its);
 else TMesh::info("No intersections detected.\n");

 FOREACHVVVERTEX(selV, v, n) UNMARK_BIT(v,5);
 if (isSelection) {delete(selT); delete(selV);}

 return its;
}

void jitterIncrease(char *f)
{
	bool isnegative = (f[0] == '-');
	int l = strlen(f);

	if (isnegative)
	{
	 for (int i = l - 1; i >= 1; i--)
		if (f[i] == '0') f[i] = '9';
		else if (f[i] == '.') continue;
		else { f[i]--; break; }
	} else
	{
	 for (int i = l - 1; i >= 0; i--)
		if (f[i] == '9') f[i] = '0';
		else if (f[i] == '.') continue;
		else { f[i]++; break; }
	}
}

void jitterDecrease(char *f)
{
	bool isnegative = (f[0] == '-');
	int l = strlen(f);

	if (isnegative)
	{
	 for (int i = l - 1; i >= 1; i--)
		if (f[i] == '9') f[i] = '0';
		else if (f[i] == '.') continue;
		else { f[i]++; break; }
	} else
	{
	 for (int i = l - 1; i >= 0; i--)
		if (f[i] == '0') f[i] = '9';
		else if (f[i] == '.') continue;
		else { f[i]--; break; }
	}
}

void jitterCoordinate(coord& c, int j)
{
 char floatver[32];
 float x;

 sprintf(floatver, "%f", TMESH_TO_FLOAT(c));
 if (j > 0) jitterIncrease(floatver);
 else if (j<0) jitterDecrease(floatver);
 sscanf(floatver, "%f", &x); c = x;
}

bool Basic_TMesh::safeCoordBackApproximation()
{
	Node *n;
	Vertex *v;

	deselectTriangles();

	FOREACHVERTEX(v, n)
	{
		jitterCoordinate(v->x, 0);
		jitterCoordinate(v->y, 0);
		jitterCoordinate(v->z, 0);
	}
	
	Edge *e;
	Vertex *ov1, *ov2;

	int pnos = 0, nos;
	nos = 0; FOREACHEDGE(e, n) if (e->overlaps()) nos++;

	do
	{
		pnos = nos;
		FOREACHEDGE(e, n) if (e->overlaps())
		{
			ov1 = e->t1->oppositeVertex(e);
			ov2 = e->t2->oppositeVertex(e);
			v = (Point::squaredTriangleArea3D(e->v1, e->v2, ov1) < Point::squaredTriangleArea3D(e->v1, e->v2, ov2)) ? (ov1) : (ov2);
			for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) for (int c = -1; c <= 1; c++)
			{
				jitterCoordinate(v->x, a); jitterCoordinate(v->y, b); jitterCoordinate(v->z, c);
				if (e->overlaps())
				{
					jitterCoordinate(v->x, -a); jitterCoordinate(v->y, -b); jitterCoordinate(v->z, -c);
				} else a = b = c = 2;
			}
		}
		nos = 0; FOREACHEDGE(e, n) if (e->overlaps()) nos++;
	} while (nos < pnos);

//	if (nos) TMesh::warning("%d overlaps could not be removed.\n", nos);
	return (nos == 0);
}








bool remints_appendCubeToList(Triangle *t0, List& l)
{
 if (!IS_VISITED(t0) || IS_BIT(t0, 6)) return false;

 Triangle *t, *s;
 Vertex *v;
 List triList(t0);
 MARK_BIT(t0,6);
 coord minx = DBL_MAX, maxx = -DBL_MAX, miny = DBL_MAX, maxy = -DBL_MAX, minz = DBL_MAX, maxz = -DBL_MAX;

 while(triList.numels())
 {
  t = (Triangle *)triList.popHead();
  v = t->v1();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  v = t->v2();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  v = t->v3();
  minx=MIN(minx,v->x); miny=MIN(miny,v->y); minz=MIN(minz,v->z);
  maxx=MAX(maxx,v->x); maxy=MAX(maxy,v->y); maxz=MAX(maxz,v->z);
  if ((s = t->t1()) != NULL && !IS_BIT(s, 6) && IS_VISITED(s)) {triList.appendHead(s); MARK_BIT(s,6);}
  if ((s = t->t2()) != NULL && !IS_BIT(s, 6) && IS_VISITED(s)) {triList.appendHead(s); MARK_BIT(s,6);}
  if ((s = t->t3()) != NULL && !IS_BIT(s, 6) && IS_VISITED(s)) {triList.appendHead(s); MARK_BIT(s,6);}
 }

 l.appendTail(new Point(minx, miny, minz));
 l.appendTail(new Point(maxx, maxy, maxz));
 return true;
}

bool remints_isVertexInCube(Vertex *v, List& loc)
{
 Node *n;
 Point *p1, *p2;
 FOREACHNODE(loc, n)
 {
  p1 = (Point *)n->data; n=n->next(); p2 = (Point *)n->data;
  if (!(v->x < p1->x || v->y < p1->y || v->z < p1->z ||
      v->x > p2->x || v->y > p2->y || v->z > p2->z)) return true;
 }

 return false;
}

void remints_selectTrianglesInCubes(Basic_TMesh *tin)
{
 Triangle *t;
 Vertex *v;
 Node *n;
 List loc;
 FOREACHVTTRIANGLE((&(tin->T)), t, n) remints_appendCubeToList(t, loc);
 FOREACHVVVERTEX((&(tin->V)), v, n) if (remints_isVertexInCube(v, loc)) MARK_BIT(v, 5);
 FOREACHVTTRIANGLE((&(tin->T)), t, n)
 {
	 UNMARK_BIT(t, 6);
	 if (IS_BIT(t->v1(), 5) || IS_BIT(t->v2(), 5) || IS_BIT(t->v3(), 5)) MARK_VISIT(t);
 }
 FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_BIT(v, 5);
 loc.freeNodes();
}


// returns true on success

bool Basic_TMesh::strongIntersectionRemoval(int max_iters)
{
 int n, iter_count = 0;
 bool qstatus = TMesh::quiet;

 TMesh::info("Removing self-intersections...\n");

 while ((++iter_count) <= max_iters && selectIntersectingTriangles())
 {
  for (n=1; n<iter_count; n++) growSelection();
  removeSelectedTriangles();
  removeSmallestComponents();
  TMesh::quiet = true; fillSmallBoundaries(E.numels(), false); TMesh::quiet = qstatus;
  coordBackApproximation();
  remints_selectTrianglesInCubes(this);
 }

 if (iter_count > max_iters) return false;
 return true;
}

} //namespace T_MESH
