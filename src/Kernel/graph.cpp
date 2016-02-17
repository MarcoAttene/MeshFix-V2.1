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

#include "graph.h"

namespace T_MESH
{

graphEdge *graphNode::getEdge(graphNode *gn)
{
 graphEdge *ge;
 Node *n = edges.head();
 while (n!=NULL)
	 if ((ge=((graphEdge *)n->data))->oppositeNode(this)==gn) return ge;
	 else n=n->next();
 return NULL;
}

graphEdge::graphEdge(graphNode *a, graphNode *b)
{
 n1=a; n2=b;
 n1->edges.appendHead(this);
 n2->edges.appendHead(this);
}


graphEdge *Graph::createEdge(graphNode *n1, graphNode *n2)
{
 Node *n;
 FOREACHNODE(n1->edges, n)
  if (((graphEdge *)n->data)->hasNode(n2))
   return (graphEdge *)n->data;

 edges.appendHead(new graphEdge(n1, n2));
 return (graphEdge *)edges.head()->data;
}


void graphEdge::collapse()
{
 Node *n;
 graphEdge *e;
 graphNode *nx;

 while ((e = (graphEdge *)n2->edges.popHead()) != NULL)
  if (e != this)
  { 
   ((e->n1 == n2)?(e->n1):(e->n2)) = n1;
   n1->edges.appendHead(e);
  }

 FOREACHNODE(n1->edges, n)
 {
  e = (graphEdge *)n->data;
  if (!e->isUnlinked()) e->oppositeNode(n1)->mask = 0;
 }

 n2->mask = 1;
 FOREACHNODE(n1->edges, n)
 {
  e = (graphEdge *)n->data;
  if (e != this)
  {
   nx = e->oppositeNode(n1);
   if (nx->mask) {nx->edges.removeNode(e); e->makeUnlinked();}
   nx->mask = 1;
  }
 }

 n = n1->edges.head();
 while (n != NULL)
 {
  e = (graphEdge *)n->data;
  n = n->next();
  if (e->isUnlinked()) n1->edges.removeCell((n!=NULL)?(n->prev()):n1->edges.tail());
 }

 FOREACHNODE(n1->edges, n)
   ((graphEdge *)n->data)->oppositeNode(n1)->mask = 0;

 n1->edges.removeNode(this);

 makeUnlinked();
}


Graph::~Graph()
{
 graphNode *gn;
 graphEdge *ge;
 while ((gn=(graphNode *)nodes.popHead())!=NULL) delete gn;
 while ((ge=(graphEdge *)edges.popHead())!=NULL) delete ge;
}

void Graph::deleteUnlinkedElements()
{
 Node *n;
 graphNode *gn;
 graphEdge *ge;

 n = nodes.head();
 while (n != NULL)
 {
  gn = (graphNode *)n->data;
  n = n->next();
  if (gn->isIsolated())
  {
   nodes.removeCell((n!=NULL)?(n->prev()):nodes.tail());
   delete(gn);
  }
 }

 n = edges.head();
 while (n != NULL)
 {
  ge = (graphEdge *)n->data;
  n = n->next();
  if (ge->isUnlinked())
  {
   edges.removeCell((n!=NULL)?(n->prev()):edges.tail());
   delete(ge);
  }
 }
}

void Graph::unlinkEdge(graphEdge *e)
{
 e->n1->edges.removeNode(e);
 e->n2->edges.removeNode(e);
 e->makeUnlinked();
}

void Graph::destroyEdge(graphEdge *e)
{
 unlinkEdge(e);
 edges.removeNode(e);
 delete e;
}

graphNode *Graph::unlinkNode(graphNode *a)
{
 graphEdge *e;
 while ((e = (graphEdge *)a->edges.popHead())!=NULL) unlinkEdge(e);
 return a;
}

void graphEdge::invert()
{
 graphNode *tmp = n1;
 n1 = n2;
 n2 = tmp;
}

bool Graph::isConnected()
{
 if (nodes.numels() < 2) return true;

 unsigned char *nmask = new unsigned char[nodes.numels()];
 Node *n;
 graphNode *p, *q;
 int i;

 for (i=0, n=nodes.head(); n!=NULL; n=n->next(), i++)
 {
  p = (graphNode *)n->data;
  nmask[i]=p->mask;
  p->mask=0;
 }

 p = (graphNode *)nodes.head()->data;
 List todo(p); p->mask = 1;
 while ((p = (graphNode *)todo.popHead())!=NULL)
 {
  for (n=p->edges.head(); n!=NULL; n=n->next())
  {
   q = ((graphEdge *)n->data)->oppositeNode(p);
   if (q->mask==0) {todo.appendTail(q); q->mask=1;}
  }
 }

 bool is_connected = true;
 for (i=0, n=nodes.head(); n!=NULL; n=n->next(), i++)
 {
  p = (graphNode *)n->data;
  if (p->mask==0) is_connected=false;
  p->mask = nmask[i];
 }

 return is_connected;
}

} //namespace T_MESH
