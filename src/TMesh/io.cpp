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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <fstream>

namespace T_MESH
{

#define VRML1_HEADER 			"#VRML V1.0 ascii"
#define VRML1_HSIZE			16
#define VRML2_HEADER 			"#VRML V2.0 utf8"
#define VRML2_HSIZE			15
#define OFF_HEADER 			"OFF"
#define OFF_HSIZE			3
#define EFF_HEADER 			"EFF"
#define EFF_HSIZE			3
#define PLY_HEADER 			"ply"
#define PLY_HSIZE			3
#define IV_HEADER 			"#Inventor V2.1 ascii"
#define IV_HSIZE			20

#define PLY_FORMAT_ASCII	0
#define PLY_FORMAT_BIN_L	1
#define PLY_FORMAT_BIN_B	2

#define TVI1(a) (TMESH_TO_INT(((Triangle *)a->data)->v1()->x))
#define TVI2(a) (TMESH_TO_INT(((Triangle *)a->data)->v2()->x))
#define TVI3(a) (TMESH_TO_INT(((Triangle *)a->data)->v3()->x))

inline void PRINT_HEADING_COMMENT(FILE *f)
{
 if (TMesh::app_name != NULL)
 {
  fprintf(f, "# File created by %s",TMesh::app_name);
  if (TMesh::app_version != NULL)
  {
   fprintf(f, " v%s",TMesh::app_version);
   if (TMesh::app_year != NULL) fprintf(f, " (%s)",TMesh::app_year);
  }
  fprintf(f, "\n");
  if (TMesh::app_url != NULL) fprintf(f, "# %s\n",TMesh::app_url);
 }
 fprintf(f, "\n");
}

inline void PRINT_PLY_COMMENT(FILE *f)
{
 if (TMesh::app_name != NULL)
 {
  fprintf(f, "comment File created by %s",TMesh::app_name);
  if (TMesh::app_version != NULL)
  {
   fprintf(f, " v%s",TMesh::app_version);
   if (TMesh::app_year != NULL) fprintf(f, " (%s)",TMesh::app_year);
  }
  fprintf(f, "\n");
  if (TMesh::app_url != NULL) fprintf(f, "comment %s\n",TMesh::app_url);
 }
}

/// Returns TRUE if the two strings are equal in a case-insensitive sense /////

inline bool sameString(const char *a, const char *b)
{
	int i = 0;
	while (a[i] != '\0' && b[i] != '\0')
	{
		if (tolower(a[i]) != tolower(b[i])) return false;
		i++;
	}

	return (a[i] == '\0' && b[i] == '\0');
}



// Swap endian-ness for four-byte elements

inline void endian_swap_long(unsigned char *p)
{
 unsigned char b0,b1,b2,b3;

 b0 = *p; b1 = *(p+1); b2 = *(p+2); b3 = *(p+3);
 *p = b3; *(p+1) = b2; *(p+2) = b1; *(p+3) = b0;
}


// Read one line (max 1024 chars) and exit if EOF

char *readLineFromFile(FILE *in, bool exit_on_eof = 1)
{
#define MAX_READLINE_CHARS	1024
 static char line[MAX_READLINE_CHARS];
 int i=0;
 char c;

 while ((c = fgetc(in)) != '\n' && i<(MAX_READLINE_CHARS-1))
   if (c==EOF)
   {
    if (exit_on_eof) TMesh::error("\nUnexpected end of file!\n");
    else return NULL;
   }
   else if (c != '\r') line[i++] = c;
 line[i] = '\0';

 if (i==MAX_READLINE_CHARS-1)
  TMesh::warning("readLineFromFile: Line is too long. Truncated !\n");

 return line;
}


// Looks for a keyword 'kw' in an ASCII file referenced through 'fp'.
// The file pointer is set to the byte right after the first keyword matched.
// Return 1 on success (keyword match), 0 otherwise.


bool seek_keyword(FILE *fp, const char *kw)
{
 static char s[256];
 s[0]='\0';
 do fscanf(fp,"%255s",s); while (strcmp(s,kw) && !feof(fp));
 if (feof(fp)) return 0;
 return 1;
}


inline void skipCommentAndBlankLines(FILE *fp)
{
 long pos0;
 char *line, s[2];
 do {pos0 = ftell(fp); line = readLineFromFile(fp);} while (line[0] == '#' || line[0] == '\0' || !sscanf(line,"%1s",s));
 fseek(fp, pos0, SEEK_SET);
}


////////////////////// Dispatch the load ///////////////////////////

int Basic_TMesh::load(const char *fname, const bool doupdate)
{
 FILE *fp;
 char header[256];
 size_t as;
 int err = IO_UNKNOWN;

 if ((fp = fopen(fname,"r")) == NULL) return IO_CANTOPEN;
 as = fread(header, 1, 256, fp);
 fclose(fp);

 if (as >= VRML1_HSIZE && !strncmp(header, VRML1_HEADER, VRML1_HSIZE)) err = loadVRML1(fname);
 else if (as >= VRML2_HSIZE && !strncmp(header, VRML2_HEADER, VRML2_HSIZE)) err = loadVRML2(fname);
 else if (as >= OFF_HSIZE && !strncmp(header, OFF_HEADER, OFF_HSIZE)) err = loadOFF(fname);
 else if (as >= EFF_HSIZE && !strncmp(header, EFF_HEADER, EFF_HSIZE)) err = loadEFF(fname);
 else if (as >= PLY_HSIZE && !strncmp(header, PLY_HEADER, PLY_HSIZE)) err = loadPLY(fname);
 else if (as >= IV_HSIZE && !strncmp(header, IV_HEADER, IV_HSIZE)) err = loadIV(fname);
 else if (sameString((char *)(fname+strlen(fname)-4), (char *)".obj")) err = loadOBJ(fname);
 else if (sameString((char *)(fname+strlen(fname)-4), (char *)".tri")) err = loadVerTri(fname);
 else if (sameString((char *)(fname+strlen(fname)-4), (char *)".stl")) err = loadSTL(fname);

 if (!err && doupdate) eulerUpdate();
 if (!err) TMesh::setFilename(fname);

 return err;
}

int Basic_TMesh::append(const char *filename, const bool doupdate)
{
 if (!T.numels()) return load(filename, doupdate);
 Basic_TMesh ntin;
 int err = ntin.load(filename, 0);

 if (err) return err;

 V.joinTailList(&(ntin.V));
 E.joinTailList(&(ntin.E));
 T.joinTailList(&(ntin.T));
 if (doupdate) eulerUpdate();
 else d_boundaries = d_handles = d_shells = 1;

 return 0;
}

int Basic_TMesh::save(const char *fname, bool back_approx)
{
 char nfname[4096];
 strcpy(nfname, fname);

 int rv;
 size_t i=strlen(fname)-1;
 while (i>0 && fname[i] != '.') i--;

 if (i==0) {strcat(nfname,".wrl"); i=strlen(fname);}

 if (sameString(nfname+i, ".wrl")) rv = saveVRML1(nfname);
 else if (sameString(nfname+i, ".iv")) rv = saveIV(nfname);
 else if (sameString(nfname + i, ".off")) rv = saveOFF(nfname);
 else if (sameString(nfname + i, ".eff")) rv = saveEFF(nfname);
 else if (sameString(nfname + i, ".ply")) rv = savePLY(nfname);
 else if (sameString(nfname+i, ".obj")) rv = saveOBJ(nfname);
 else if (sameString(nfname+i, ".stl")) rv = saveSTL(nfname);
 else if (sameString(nfname+i, ".tri")) {nfname[i]='\0'; rv = saveVerTri(nfname);}
 else
 {
  TMesh::warning("Unknown extension '%s'.\n",nfname+i);
  TMesh::warning("I did not save anything.\n");
  TMesh::warning("Recognized extensions are:");
  TMesh::warning(".wrl (ASCII VRML 1.0)\n");
  TMesh::warning(".iv (Open Inventor 2.1)\n");
  TMesh::warning(".off (Object File Format)\n");
  TMesh::warning(".eff (Exact File Format)\n");
  TMesh::warning(".obj (Wavefront/Java3D)\n");
  TMesh::warning(".stl (Stereolithography)\n");
  TMesh::warning(".ply (Ascii PLY 1.0 Format)\n");
  TMesh::warning(".tri (IMATI ver-tri File Format)\n");
  return 0;
 }

 if (!rv && back_approx) coordBackApproximation();

 return rv;
}


// This part is common to all the loaders

Triangle * Basic_TMesh::CreateTriangleFromVertices(ExtVertex *vari1, ExtVertex *vari2, ExtVertex *vari3)
{
 Edge *e1, *e2, *e3;
 Triangle *t = NULL;

 e1 = CreateEdge(vari1,vari2); if (e1->t1 != NULL && e1->t2 != NULL) MARK_BIT(e1,5);
 e2 = CreateEdge(vari2,vari3); if (e2->t1 != NULL && e2->t2 != NULL) MARK_BIT(e2,5);
 e3 = CreateEdge(vari3,vari1); if (e3->t1 != NULL && e3->t2 != NULL) MARK_BIT(e3,5);
 if (IS_BIT(e1,5)) {e1 = CreateEdge(vari1,vari2,0); MARK_BIT(e1,5);}
 if (IS_BIT(e2,5)) {e2 = CreateEdge(vari2,vari3,0); MARK_BIT(e2,5);}
 if (IS_BIT(e3,5)) {e3 = CreateEdge(vari3,vari1,0); MARK_BIT(e3,5);}

 if ((t=CreateUnorientedTriangle(e1,e2,e3)) == NULL)
 {
  if (e3->t1 == NULL && e3->t2 == NULL)
  {
   E.freeNode(e3);
   vari3->VE.removeNode(e3); vari1->VE.removeNode(e3);
   if (vari3->v->e0 == e3) vari3->v->e0 = NULL;
   if (vari1->v->e0 == e3) vari1->v->e0 = NULL;
  }
  if (e2->t1 == NULL && e2->t2 == NULL)
  {
   E.freeNode(e2);
   vari2->VE.removeNode(e2); vari3->VE.removeNode(e2);
   if (vari2->v->e0 == e2) vari2->v->e0 = NULL;
   if (vari3->v->e0 == e2) vari3->v->e0 = NULL;
  }
  if (e1->t1 == NULL && e1->t2 == NULL)
  {
   E.freeNode(e1);
   vari1->VE.removeNode(e1); vari2->VE.removeNode(e1);
   if (vari1->v->e0 == e1) vari1->v->e0 = NULL;
   if (vari2->v->e0 == e1) vari2->v->e0 = NULL;
  }
 }

 return t;
}


// This part is common to all the loaders

Triangle * Basic_TMesh::CreateIndexedTriangle(ExtVertex **var, int i1, int i2, int i3)
{
 return CreateTriangleFromVertices(var[i1], var[i2], var[i3]);
}


// This part is common to all the loaders

void Basic_TMesh::closeLoadingSession(FILE *fp, int loaded_faces, ExtVertex **var, bool triangulate)
{
 int i, nv = V.numels();

 fclose(fp);

 if (var != NULL)
 {
  for (i=0; i<nv; i++) delete(var[i]);
  free(var);
 }

 if (loaded_faces)
 {
  TMesh::info("Loaded %d vertices and %d faces.\n",nv,loaded_faces);
  if (triangulate) TMesh::warning("Some polygonal faces needed to be triangulated.\n");
  fixConnectivity();
 }

 d_boundaries = d_handles = d_shells = 1;
}


// This method should be called after a Save to ascii file to ensure
// coherence between in-memory data and saved data.

void Basic_TMesh::coordBackApproximation()
{
 Node *n;
 Vertex *v;
 char floatver[32];
 float x;

 FOREACHVERTEX(v, n)
 {
  sprintf(floatver, "%f", TMESH_TO_FLOAT(v->x)); sscanf(floatver, "%f", &x); v->x = x;
  sprintf(floatver, "%f", TMESH_TO_FLOAT(v->y)); sscanf(floatver, "%f", &x); v->y = x;
  sprintf(floatver, "%f", TMESH_TO_FLOAT(v->z)); sscanf(floatver, "%f", &x); v->z = x;
 }
}

////////////////////// Loads VRML 1.0 format ///////////////////////////

int Basic_TMesh::loadVRML1(const char *fname)
{
 FILE *fp;
 Node *n;
 float x,y,z;
 int i,i1,i2,i3,i4,nv=0,triangulate=0;
 Vertex *v;

 coord cx;

 if ((fp = fopen(fname,"r")) == NULL) return IO_CANTOPEN;

 if (!seek_keyword(fp, "point")) {closeLoadingSession(fp, 0, NULL, 0); return IO_FORMAT;}
 if (!seek_keyword(fp, "[")) {closeLoadingSession(fp, 0, NULL, 0); return IO_FORMAT;}

 while (fscanf(fp, "%f %f %f,", &x, &y, &z) == 3) { V.appendTail(newVertex(x, y, z)); cx = coord(x); /*printf("%f %f\n", x, cx.toFloat());*/ }
 nv = V.numels();
 ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);
 i = 0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

 if (!seek_keyword(fp, "coordIndex")) {closeLoadingSession(fp, 0, var, 0); return IO_FORMAT;}
 if (!seek_keyword(fp, "[")) {closeLoadingSession(fp, 0, var, 0); return IO_FORMAT;}

 i=0; TMesh::begin_progress();
 while (fscanf(fp,"%d, %d, %d,",&i1,&i2,&i3) == 3)
 {
  if (((i++)%1000) == 0) TMesh::report_progress("Loading ..%d%%",(i*100)/(nv*2));
  if (i1<0 || i2<0 || i3<0 || i1>(nv-1) || i2>(nv-1) || i3>(nv-1))
	  TMesh::error("\nloadVRML1: Invalid indices %d %d %d!\n",i1,i2,i3);
  do
  {
   if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadVRML1: Coincident indexes at face %d! Skipping.\n",i);
   else if (!CreateIndexedTriangle(var, i1, i2, i3)) TMesh::warning("\nloadVRML1: This shouldn't happen!!! Skipping triangle.\n");
   if (fscanf(fp,"%d,",&i4) != 1) TMesh::error("loadVRML1: Unexpected end of file at face %d!\n",i);
   i2=i3; i3=i4;
   if (i4 != -1) triangulate=1;
  } while (i4 != -1);
 }
 TMesh::end_progress();

 closeLoadingSession(fp, i, var, (triangulate != 0));
 TMesh::setFilename(fname);

 return 0;
}


int Basic_TMesh::loadIV(const char *fname)
{
 return loadVRML1(fname);
}


////////////////////// Loads OFF format ///////////////////////////

int Basic_TMesh::loadOFF(const char *fname)
{
 FILE *fp;
 Node *n;
 char s[256], *line;
 float x,y,z;
 int i,j,i1,i2,i3,i4,nv,nt,ne,triangulate=0;
 Vertex *v;

 if ((fp = fopen(fname,"rb")) == NULL) return IO_CANTOPEN;

 fscanf(fp,"%255s",s);
 if (strcmp(s,"OFF") || feof(fp)) return IO_FORMAT;
 do {line = readLineFromFile(fp);} while (line[0] == '#' || line[0] == '\0' || !sscanf(line,"%256s",s));
 if (sscanf(line,"%d %d %d",&nv,&nt,&ne) < 3) return IO_FORMAT;
 if (nv < 3) TMesh::error("\nloadOFF: Sorry. Can't load objects with less than 3 vertices.\n");
 if (nt < 1) TMesh::error("\nloadOFF: Sorry. Can't load objects with no faces.\n");

 skipCommentAndBlankLines(fp);

 for (i = 0; i<nv; i++)
  if (fscanf(fp,"%f %f %f",&x,&y,&z) == 3) V.appendTail(newVertex(x,y,z));
  else TMesh::error("\nloadOFF: Couldn't read coordinates for vertex # %d\n",i);

 ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);
 i=0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

 skipCommentAndBlankLines(fp);

 TMesh::begin_progress();
 for (i=0; i<nt; i++)
 {
  if (fscanf(fp,"%d %d %d %d",&i4,&i1,&i2,&i3) == 4)
  {
   if ((i%1000) == 0) TMesh::report_progress("Loading ..%d%%",(i*100)/(nv*2));
   if (i1<0 || i2<0 || i3<0 || i4<3 || i1>(nv-1) || i2>(nv-1) || i3>(nv-1)) TMesh::error("\nloadOFF: Invalid index at face %d!\n",i);
   for (j=3; j<=i4; j++)
   {
    if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadOFF: Coincident indexes at triangle %d! Skipping.\n",i);
    else if (!CreateIndexedTriangle(var, i1, i2, i3)) TMesh::warning("\nloadOFF: This shouldn't happen!!! Skipping triangle.\n");
    i2 = i3;
    if (j<i4)
    {
     if (fscanf(fp,"%d",&i3) != 1) TMesh::error("\nloadOFF: Couldn't read indexes for face # %d\n",i);
     else triangulate=1;
    }
   }
  }
  else TMesh::error("\nloadOFF: Couldn't read indexes for face # %d\n",i);
 }

 TMesh::end_progress();

 closeLoadingSession(fp, i, var, (triangulate != 0));
 TMesh::setFilename(fname);

 return 0;
}

////////////////////// Loads EFF format ///////////////////////////

int Basic_TMesh::loadEFF(const char *fname)
{
	std::ifstream is;
	is.open(fname);
	if (!is.is_open()) return IO_CANTOPEN;
	Node *n;
	char s[256];
	coord x, y, z;
	int i, i1, i2, i3, nv=-1, nt=-1, triangulate = 0;
	Vertex *v;


	is >> s;
	if (is.eof() || is.fail() || strcmp(s, "EFF")) return IO_FORMAT;
	is >> nv; if (is.eof() || is.fail()) return IO_FORMAT;
	is >> nt; if (is.eof() || is.fail()) return IO_FORMAT;
	if (nv < 3) TMesh::error("\nloadOFF: Sorry. Can't load objects with less than 3 vertices.\n");
	if (nt < 1) TMesh::error("\nloadOFF: Sorry. Can't load objects with no faces.\n");

	TMesh::useRationals(true);
	for (i = 0; i < nv; i++)
	{
		is >> x >> y >> z;
		if (is.eof() || is.fail()) TMesh::error("\nloadEFF: Couldn't read coordinates for vertex # %d\n", i);
		V.appendTail(newVertex(x, y, z));
	}

	ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);
	i = 0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

	TMesh::begin_progress();
	for (i = 0; i<nt; i++)
	{
		is >> i1 >> i2 >> i3;
		if (is.eof() || is.fail()) TMesh::error("\nloadEFF: Couldn't read indexes for face # %d\n", i);
		else
		{
			if ((i % 1000) == 0) TMesh::report_progress("Loading ..%d%%", (i * 100) / (nv * 2));
			if (i1<0 || i2<0 || i3<0 || i1>(nv - 1) || i2>(nv - 1) || i3>(nv - 1)) TMesh::error("\nloadEFF: Invalid index at face %d!\n", i);
			if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadEFF: Coincident indexes at triangle %d! Skipping.\n", i);
			else if (!CreateIndexedTriangle(var, i1, i2, i3)) TMesh::warning("\nloadEFF: This shouldn't happen!!! Skipping triangle.\n");
		}
	}

	TMesh::end_progress();

	is.close();

	for (i = 0; i<nv; i++) delete(var[i]);
	free(var);

	TMesh::info("Loaded %d vertices and %d faces.\n", nv, nt);
	fixConnectivity();

	d_boundaries = d_handles = d_shells = 1;
	TMesh::setFilename(fname);

	return 0;
}

////////////////////// Loads VRML 2.0 format ///////////////////////////

int Basic_TMesh::loadVRML2(const char *fname)
{
 FILE *fp;
 Node *n;
 float x,y,z;
 int i,i1,i2,i3,i4,nv=0,triangulate=0;
 Vertex *v;

 if ((fp = fopen(fname,"r")) == NULL) return IO_CANTOPEN;

 if (!seek_keyword(fp, "point")) {fclose(fp); return IO_FORMAT;}
 if (!seek_keyword(fp, "[")) {fclose(fp); return IO_FORMAT;}

 while (fscanf(fp,"%f %f %f,",&x,&y,&z) == 3) V.appendTail(newVertex(x,y,z));
 nv = V.numels();

 if (!seek_keyword(fp, "coordIndex")) {fclose(fp); return IO_FORMAT;}
 if (!seek_keyword(fp, "[")) {fclose(fp); return IO_FORMAT;}

 ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);
 i=0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

 i=0;
 TMesh::begin_progress();
 while (fscanf(fp,"%d, %d, %d,",&i1,&i2,&i3) == 3)
 {
  if (((i++)%1000) == 0) TMesh::report_progress("Loading ..%d%%",(i*100)/(nv*2));
  if (i1<0 || i2<0 || i3<0 || i1>(nv-1) || i2>(nv-1) || i3>(nv-1)) TMesh::error("\nloadVRML2: Invalid index at face %d!\n",i);
  do
  {
   if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadVRML2: Coincident indexes at triangle %d! Skipping.\n",i);
   else if (!CreateIndexedTriangle(var, i1, i2, i3)) TMesh::warning("\nloadVRML2: This shouldn't happen!!! Skipping triangle.\n");
   if (fscanf(fp,"%d,",&i4) != 1) TMesh::error("loadVRML2: Unexpected end of file at triangle %d!\n",i);
   i2=i3; i3=i4;
   if (i4 != -1) triangulate=1;
  } while (i4 != -1);
 }
 TMesh::end_progress();

 closeLoadingSession(fp, i, var, (triangulate != 0));
 TMesh::setFilename(fname);

 return 0;
}


////////////////////// Loads Ver-Tri format ////////////////////

int Basic_TMesh::loadVerTri(const char *fname)
{
 FILE *fpv, *fpt;
 int numvers, numtris, i, i1, i2, i3, a1, a2, a3;
 float x,y,z;
 char vername[256], triname[256];
 Node *n;
 Vertex *v;

 if (!sameString((char *)(fname+strlen(fname)-4), (char *)".tri")) return IO_UNKNOWN;

 strcpy(triname,fname);
 strcpy(vername,fname); vername[strlen(vername)-4]='\0';
 strcat(vername,".ver");

 if ((fpv = fopen(vername,"r")) == NULL)
 {
  fprintf(stderr,"Can't open '%s' for input !\n",vername);
  return 1;
 }
 if ((fpt = fopen(triname,"r")) == NULL)
 {
  fclose(fpv);
  fprintf(stderr,"Can't open '%s' for input !\n",triname);
  return 1;
 }

 if (!fscanf(fpv,"%d\n",&numvers) || numvers < 3) {fclose(fpv); fclose(fpt); return IO_FORMAT;}
 if (!fscanf(fpt,"%d\n",&numtris) || numtris < 1) {fclose(fpv); fclose(fpt); return IO_FORMAT;}

 for (i=0; i<numvers; i++)
  if (fscanf(fpv,"%f %f %f\n",&x,&y,&z) != 3) TMesh::error("Couldn't read %d'th vertex!\n",i+1);
  else V.appendTail(newVertex(x,y,z));
 fclose(fpv);

 ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*numvers);
 i=0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

 TMesh::begin_progress();
 for (i=0; i<numtris; i++)
  if (fscanf(fpt,"%d %d %d %d %d %d",&i1,&i2,&i3,&a1,&a2,&a3) == 6)
  {
   if (((i)%1000) == 0) TMesh::report_progress("Loading ..%d%%",(i*100)/(numtris));
   if (i1 < 1 || i2 < 1 || i3 < 1) TMesh::error("\nloadVerTri: Illegal index at triangle %d!\n",i);
   else if (i1 > (numvers) || i2 > (numvers) || i3 > (numvers)) TMesh::error("\nloadVerTri: Index out of bounds at triangle %d!\n",i);
   else if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadVerTri: Coincident indexes at triangle %d! Skipping.\n",i);
   else if (!CreateIndexedTriangle(var, i1-1, i2-1, i3-1)) TMesh::warning("\nloadVerTri: This shouldn't happen!!! Skipping triangle.\n");
  }
  else TMesh::error("loadVerTri: Couldn't read %dth triangle !\n",i+1);

 TMesh::end_progress();

 closeLoadingSession(fpt, T.numels(), var, 0);
 TMesh::setFilename(fname);

 return 0;
}


////////////////////// Saves IV 2.1 format ////////////////////

int Basic_TMesh::saveIV(const char *fname)
{
 FILE *fp;
 int i;
 char triname[256];
 Node *n;
 coord *ocds;
 Vertex *v;

 strcpy(triname,fname);

 if ((fp = fopen(triname,"w")) == NULL)
 {
  TMesh::warning("Can't open '%s' for output !\n",triname);
  return 1;
 }

 fprintf(fp,"#Inventor V2.1 ascii\n\n");
 PRINT_HEADING_COMMENT(fp);
 fprintf(fp,"Separator {\n");
 fprintf(fp," Coordinate3 {\n  point [\n");

 FOREACHVERTEX(v, n) fprintf(fp, "   %f %f %f,\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));

 fprintf(fp,"  ]\n }\n");
 fprintf(fp," IndexedFaceSet {\n  coordIndex [\n");

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
 i=0; FOREACHVERTEX(v, n) v->x = i++;

 FOREACHNODE(T, n) fprintf(fp,"   %d, %d, %d, -1,\n",TVI1(n),TVI2(n),TVI3(n));

 fprintf(fp,"  ]\n }\n");
 fprintf(fp,"}\n");

 fclose(fp);
 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}


////////////////////// Saves VRML 1.0 format ////////////////////

int Basic_TMesh::saveVRML1(const char *fname, const int mode)
{
 FILE *fp;
 int i;
 unsigned int pkc;
 char triname[256];
 Node *n;
 Vertex *v;
 Triangle *t;
 coord *ocds;

 strcpy(triname,fname);

 if ((fp = fopen(triname,"w")) == NULL)
 {
  TMesh::warning("Can't open '%s' for output !\n",triname);
  return 1;
 }

 fprintf(fp,"#VRML V1.0 ascii\n\n");
 PRINT_HEADING_COMMENT(fp);
 fprintf(fp,"Separator {\n");
 fprintf(fp," Coordinate3 {\n  point [\n");

 FOREACHVERTEX(v, n) fprintf(fp, "   %f %f %f,\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));

 fprintf(fp,"  ]\n }\n");

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) {ocds[i] = v->x; v->x = i++;}

 switch (mode)
 {
  case IO_CSAVE_OVERALL:
   fprintf(fp,"Material {\n diffuseColor 0.6 0.6 0.6\n}\n");
   break;
  case IO_CSAVE_PERFACE:
   fprintf(fp,"Material {\n diffuseColor [\n");
   FOREACHTRIANGLE(t, n)
   {
    pkc = (unsigned int)((j_voidint)t->info);
    fprintf(fp,"  %f %f %f,\n",((pkc>>24)&0x000000ff)/255.0,((pkc>>16)&0x000000ff)/255.0,((pkc>>8)&0x000000ff)/255.0);
   }
   fprintf(fp," ]\n}\nMaterialBinding {\n value PER_FACE_INDEXED\n}\n");
   break;
  case IO_CSAVE_PERVERTEX:
   fprintf(fp,"Material {\n diffuseColor [\n");
   FOREACHVERTEX(v, n)
   {
    pkc = (unsigned int)((j_voidint)v->info);
    fprintf(fp,"  %f %f %f,\n",((pkc>>24)&0x000000ff)/255.0,((pkc>>16)&0x000000ff)/255.0,((pkc>>8)&0x000000ff)/255.0);
   }
   fprintf(fp," ]\n}\nMaterialBinding {\n value PER_VERTEX_INDEXED\n}\n");
   break;
  case IO_CSAVE_PERFACE_INDEXED:
   fprintf(fp,"Material {\n diffuseColor [\n");
   fprintf(fp,"1.0 1.0 1.0,\n1.0 0.0 0.0,\n0.0 1.0 0.0,\n0.0 0.0 1.0,\n 0.8 0.8 0.0\n");
   fprintf(fp," ]\n}\nMaterialBinding {\n value PER_FACE_INDEXED\n}\n");
   break;
  case IO_CSAVE_PERVERTEX_INDEXED:
   fprintf(fp,"Material {\n diffuseColor [\n");
   fprintf(fp,"1.0 1.0 1.0,\n1.0 0.0 0.0,\n0.0 1.0 0.0,\n0.0 0.0 1.0,\n 0.8 0.8 0.0\n");
   fprintf(fp," ]\n}\nMaterialBinding {\n value PER_VERTEX_INDEXED\n}\n");
   break;
  default: TMesh::error("Basic_TMesh::saveVRML1. Unknown mode %d\n",mode);
 }

 fprintf(fp," IndexedFaceSet {\n  coordIndex [\n");

 FOREACHTRIANGLE(t, n)
	 fprintf(fp, "   %d, %d, %d, -1,\n", TMESH_TO_INT(t->v1()->x), TMESH_TO_INT(t->v2()->x), TMESH_TO_INT(t->v3()->x));

 fprintf(fp,"  ]\n");

 if (mode != IO_CSAVE_OVERALL)
 {
  fprintf(fp,"  materialIndex [\n");
  switch (mode)
  {
   case IO_CSAVE_PERFACE_INDEXED:
    FOREACHTRIANGLE(t, n) fprintf(fp,"   %d,\n",t->mask);
    break;
   case IO_CSAVE_PERVERTEX_INDEXED:
    FOREACHTRIANGLE(t, n) fprintf(fp,"   %d, %d, %d, -1,\n",t->v1()->mask,t->v2()->mask,t->v3()->mask);
    break;
   case IO_CSAVE_PERFACE:
    i=0; FOREACHTRIANGLE(t, n) fprintf(fp,"   %d,\n",i++);
    break;
   case IO_CSAVE_PERVERTEX:
   FOREACHTRIANGLE(t, n) fprintf(fp, "   %d, %d, %d, -1,\n", TMESH_TO_INT(t->v1()->x), TMESH_TO_INT(t->v2()->x), TMESH_TO_INT(t->v3()->x));
	break;
  }
  fprintf(fp,"  ]\n");
 }

 fprintf(fp," }\n}\n");

 fclose(fp);
 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}


////////////////////// Saves OFF format ///////////////////////////

int Basic_TMesh::saveOFF(const char *fname)
{
 FILE *fp;
 int i;
 char triname[256];
 Node *n;
 coord *ocds;
 Vertex *v;

 strcpy(triname,fname);

 if ((fp = fopen(triname,"w")) == NULL)
 {
  TMesh::warning("Can't open '%s' for output !\n",triname);
  return 1;
 }

 fprintf(fp,"OFF\n");
 PRINT_HEADING_COMMENT(fp);
 fprintf(fp,"%d %d 0\n",V.numels(),T.numels());

 FOREACHVERTEX(v, n) fprintf(fp, "%f %f %f\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
 i=0; FOREACHVERTEX(v, n) v->x = i++;

 FOREACHNODE(T, n) fprintf(fp,"3 %d %d %d\n",TVI1(n),TVI2(n),TVI3(n));

 fclose(fp);
 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}


////////////////////// Saves Ver-Tri format ////////////////////

//#define SAVE_INFO
int Basic_TMesh::saveVerTri(const char *fname)
{
#ifdef SAVE_INFO
 TMesh::warning("saveVerTri: Assuming that the vertex info field is allocated!\n");
 FILE *fpj;
 char jkkname[256];
#endif
 FILE *fpv, *fpt;
 int i, i1, i2, i3, a1, a2, a3;
 char vername[256], triname[256];
 Node *n;
 Vertex *v;
 Triangle *t, *t1, *t2, *t3;
 coord *ocds;

 strcpy(triname,fname);
 strcpy(vername,fname);
 strcat(triname,".tri");
 strcat(vername,".ver");

#ifdef SAVE_INFO
 strcpy(jkkname,fname);
 strcat(jkkname,".jkk");
#endif

 if ((fpv = fopen(vername,"w")) == NULL)
 {
  fprintf(stderr,"Can't open '%s' for output !\n",vername);
  return 1;
 }
 if ((fpt = fopen(triname,"w")) == NULL)
 {
  fclose(fpv);
  fprintf(stderr,"Can't open '%s' for output !\n",triname);
  return 1;
 }
#ifdef SAVE_INFO
 if ((fpj = fopen(jkkname,"w")) == NULL)
 {
  fclose(fpv); fclose(fpt);
  fprintf(stderr,"Can't open '%s' for output !\n",jkkname);
  return 1;
 }
#endif

 fprintf(fpv,"%d\n",V.numels());
 FOREACHVERTEX(v, n)
 {
	 fprintf(fpv, "%f %f %f\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));
 }
 fclose(fpv);

#ifdef SAVE_INFO
 for (n=V.tail; n != NULL; n=n->prev)
 {
  v = ((Vertex *)n->data);
  fprintf(fpj,"%f\n",(*((double *)(v->info))));
 }
 fclose(fpj);
#endif

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
 i=0; FOREACHVERTEX(v, n) v->x = ++i;
 i=0; FOREACHTRIANGLE(t, n) {i++; t->info = (void *)(intptr_t)i;}

 fprintf(fpt,"%d\n",T.numels());
 FOREACHTRIANGLE(t, n)
 {
  i1 = TMESH_TO_INT(t->v1()->x); i2 = TMESH_TO_INT(t->v2()->x); i3 = TMESH_TO_INT(t->v3()->x);
  t1 = t->t1(); t2 = t->t2(); t3 = t->t3();
  a1 = (t1)?((long int)(t1->info)):(0); a2 = (t2)?((long int)(t2->info)):(0); a3 = (t3)?((long int)(t3->info)):(0);
  fprintf(fpt,"%d %d %d %d %d %d\n",i1, i2, i3, a1, a2, a3);
 }
 fclose(fpt);

 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}


// Implements the cutting and stitching procedure to convert to manifold mesh //
// Assumes that singular edges to be cut and stitched are marked as BIT5.  //

bool Basic_TMesh::pinch(Edge *e1, bool with_common_vertex)
{
	List *ee = (List *)e1->info;
	if (ee == NULL) return false;
	Node *n = NULL;
	Edge *e2=NULL;
	List *ve;

	if (with_common_vertex)
	{
		e1->v1->e0 = e1; ve = e1->v1->VE();
		FOREACHVEEDGE(ve, e2, n) if (e2 != e1 && e2->isOnBoundary() && (*(e2->oppositeVertex(e1->v1))) == (*(e1->v2)) && e1->merge(e2)) break;
		delete ve;
		if (n == NULL)
		{
			e1->v2->e0 = e1; ve = e1->v2->VE();
			FOREACHVEEDGE(ve, e2, n) if (e2 != e1 && e2->isOnBoundary() && (*(e2->oppositeVertex(e1->v2))) == (*(e1->v1)) && e1->merge(e2)) break;
			delete ve;
		}
	}
	else //if (ee->numels()==2)
	{
		if (e1->t1 != NULL)
		{
			FOREACHVEEDGE(ee, e2, n) if (e2 != e1 && (((*(e2->v1)) == (*(e1->v1)) && e2->t2 != NULL) || ((*(e2->v1)) == (*(e1->v2)) && e2->t1 != NULL)) && e1->merge(e2)) break;
		}
		else
		{
			FOREACHVEEDGE(ee, e2, n) if (e2 != e1 && (((*(e2->v1)) == (*(e1->v1)) && e2->t1 != NULL) || ((*(e2->v1)) == (*(e1->v2)) && e2->t2 != NULL)) && e1->merge(e2)) break;
		}
	}
	if (n == NULL) return false;

	ee->removeNode(e1); ee->removeNode(e2); e1->info = e2->info = NULL;
	if (ee->numels() == 0) delete ee;

	Edge *e, *e_1 = NULL, *e_2 = NULL;
	ve = e1->v1->VE();
	for (n = ve->head(); n != NULL; n = n->next()) if ((e = (Edge *)n->data)->info != NULL) { e_1 = e; break; }
	for (n = ve->tail(); n != NULL; n = n->prev()) if ((e = (Edge *)n->data)->info != NULL)
	{
		if ((*(e->oppositeVertex(e1->v1))) != (*(e_1->oppositeVertex(e1->v1)))) e_1 = NULL;
		break;
	}
	delete ve;

	ve = e1->v2->VE();
	for (n = ve->head(); n != NULL; n = n->next()) if ((e = (Edge *)n->data)->info != NULL) { e_2 = e; break; }
	for (n = ve->tail(); n != NULL; n = n->prev()) if ((e = (Edge *)n->data)->info != NULL)
	{
		if ((*(e->oppositeVertex(e1->v2))) != (*(e_2->oppositeVertex(e1->v2)))) e_2 = NULL;
		break;
	}
	delete ve;

	if (e_1 != NULL) pinch(e_1, true);
	if (e_2 != NULL) pinch(e_2, true);

	return true;
}

Edge *Basic_TMesh::duplicateEdge(Edge *e1)
{
	if (e1->t1 == NULL || e1->t2 == NULL) return NULL;
	Edge *e2 = newEdge(e1); //e2->invert();
	E.appendHead(e2);
	e1->t2->replaceEdge(e1, e2);
	e2->t2 = e1->t2; e1->t2 = NULL;
	return e2;
}

int Basic_TMesh::cutAndStitch()
{
	Edge *e1, *e2;
	Node *n;
	List singular_edges;

	FOREACHEDGE(e1, n) if (IS_BIT(e1, 5) && ((e2 = duplicateEdge(e1)) != NULL)) MARK_BIT(e2, 5);

	FOREACHEDGE(e1, n) if (IS_BIT(e1, 5))
	{
		singular_edges.appendHead(e1);
		UNMARK_BIT(e1, 5);
	}

	forceNormalConsistence();
	duplicateNonManifoldVertices();

	singular_edges.sort(&lexEdgeCompare);
	FOREACHEDGE(e1, n) e1->info = NULL;
	e2 = NULL;
	FOREACHVEEDGE((&singular_edges), e1, n)
	{
		if (e2 == NULL || lexEdgeCompare(e1, e2) != 0) { e1->info = new List(); e2 = e1; }
		((List *)e2->info)->appendTail(e1);
		e1->info = e2->info;
	}
	// Now each edge is either 'regular' or has the info field pointing to a list of coincident boundary edges

	// First, pinch bounded chains of singular edges starting from one endpoint
	FOREACHVEEDGE((&singular_edges), e1, n) if (e1->isLinked()) pinch(e1, true);

	// Then, pinch the remaining unbounded chains starting from any of the edges
	FOREACHVEEDGE((&singular_edges), e1, n) if (e1->isLinked()) pinch(e1, false);

	removeUnlinkedElements();

	d_boundaries = d_handles = d_shells = 1;

	return singular_edges.numels();
}


//int Basic_TMesh::cutAndStitch()
//{
//	Edge *e1, *e2;
//	Node *n;
//	List cut;
//	int i;
//
//	FOREACHEDGE(e1, n) if (IS_BIT(e1, 5))
//	{
//		if (e1->t1 != NULL && e1->t2 != NULL)
//		{
//			e2 = newEdge(e1);
//			E.appendHead(e2);
//			e1->t2->replaceEdge(e1, e2);
//			e2->t2 = e1->t2; e1->t2 = NULL;
//		}
//		cut.appendHead(e1);
//		UNMARK_BIT(e1, 5);
//	}
//
//	do
//	{
//		i = 0;
//		FOREACHVEEDGE((&cut), e1, n) if (e1->v1 != NULL) i += e1->stitch();
//	} while (i);
//
//	removeEdges();
//
//	d_boundaries = d_handles = d_shells = 1;
//
//	return cut.numels();
//}


////////////////////// PLY LOADER //////////////////////////////////////////////

int ply_parseElements(FILE *in, const char *elname)
{
 char c, keyword[64];
 int num;
 // skip comments
 if (!fscanf(in,"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");
 while (!strcmp(keyword,"comment") || !strcmp(keyword,"obj_info"))
 {
  while ((c = fgetc(in)) != '\n') if (c==EOF) TMesh::error("\nUnexpected end of file!\n");
  if (!fscanf(in,"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");
 }
 if (strcmp(keyword,"element")) TMesh::error("element definition expected!\n");
 if (!fscanf(in,"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");
 if (strcmp(keyword,elname)) TMesh::error("Sorry. Element type '%s' is not supported!\n",keyword);
 if (!fscanf(in,"%d\n",&num)) TMesh::error("Unexpected token or end of file!\n");
 if (num <= 0) TMesh::error("Unexpected empty element list!\n");

 return num;
}

void ply_checkVertexProperties(FILE *in)
{
 char keyword[64], dtype[64], dval[64];
 if (fscanf(in,"%64s %64s %64s\n",keyword,dtype,dval) < 3) TMesh::error("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) TMesh::error("property definition expected!\n");
 if (strcmp(dtype,"float") && strcmp(dtype,"float32")) TMesh::error("float property expected!\n");
 if (strcmp(dval,"x")) TMesh::error("'x' float property expected!\n");
 if (fscanf(in,"%64s %64s %64s\n",keyword,dtype,dval) < 3) TMesh::error("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) TMesh::error("property definition expected!\n");
 if (strcmp(dtype,"float") && strcmp(dtype,"float32")) TMesh::error("float property expected!\n");
 if (strcmp(dval,"y")) TMesh::error("'y' float property expected!\n");
 if (fscanf(in,"%64s %64s %64s\n",keyword,dtype,dval) < 3) TMesh::error("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) TMesh::error("property definition expected!\n");
 if (strcmp(dtype,"float") && strcmp(dtype,"float32")) TMesh::error("float property expected!\n");
 if (strcmp(dval,"z")) TMesh::error("'z' float property expected!\n");
}

int ply_getOverhead(FILE *in, int format, const char *element)
{
 char keyword[64], ptype[64], pname[64];
 int oh = 0;
 long pos = ftell(in);
 char *rline = readLineFromFile(in);
 if (!sscanf(rline,"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");
 while (!strcmp(keyword, "property"))
 {
  if (sscanf(rline,"%64s %64s %64s",keyword,ptype,pname) < 3) TMesh::error("Unexpected token or end of file!\n");
  if (!strcmp(element,"vertex") && !strcmp(pname,"x")) break;
  else if (!strcmp(element,"face") && !strcmp(ptype,"list")) break;
  pos = ftell(in);
  if (!strcmp(ptype, "char") || !strcmp(ptype, "uchar")) oh += (format)?(1):1;
  else if (!strcmp(ptype, "short") || !strcmp(ptype, "ushort")) oh += (format)?(2):1;
  else if (!strcmp(ptype, "int") || !strcmp(ptype, "uint") ||
           !strcmp(ptype, "float") || !strcmp(ptype,"float32")) oh += (format)?(4):1;
  else if (!strcmp(ptype, "double")) oh += (format)?(8):1;
  else if (!strcmp(ptype, "list")) TMesh::error("list properties other than face indices are not supported!\n");
  else TMesh::error("Unrecognized property type!\n");
  if (!sscanf(readLineFromFile(in),"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");
 }
 fseek(in, pos, SEEK_SET);

 return oh;
}

void ply_checkFaceProperties(FILE *in)
{
 char keyword[64], ltype[64], uctype[64], dtype[64], dval[64];
 if (fscanf(in,"%64s %64s %64s %64s %64s\n",keyword,ltype,uctype,dtype,dval) < 5) TMesh::error("Unexpected token or end of file!\n");
 if (strcmp(keyword,"property")) TMesh::error("property definition expected!\n");
 if (strcmp(ltype,"list")) TMesh::error("list property expected!\n");
 if (strcmp(uctype,"uchar") && strcmp(uctype,"uint8")) TMesh::error("uchar property expected!\n");
 if (strcmp(dtype,"int") && strcmp(dtype,"int32")) TMesh::error("int property expected!\n");
 if (strcmp(dval,"vertex_indices")) TMesh::error("vertex_indices property expected!\n");
}

void ply_readOverhead(FILE *in, int format, int oh)
{
 int i;
 static char token[1024];
 if (format == PLY_FORMAT_ASCII) for (i=0; i<oh; i++) fscanf(in, "%s", token);
 else for (i=0; i<oh; i++) fgetc(in);
}


int ply_readVCoords(FILE *in, int format, int ph, int oh, float *x, float *y, float *z)
{
 float vc[3];

 ply_readOverhead(in, format, ph);

 if (format == PLY_FORMAT_ASCII)
 {
  if (fscanf(in,"%f %f %f", x, y, z) < 3) TMesh::error("Unexpected token or end of file!\n");
 }
 else
 {
  if (fread(vc, 4, 3, in) < 3) TMesh::error("Unexpected end of file!\n");
  *x = vc[0]; *y = vc[1]; *z = vc[2];

  if (format == PLY_FORMAT_BIN_B)
  {
   endian_swap_long((unsigned char *)(x));
   endian_swap_long((unsigned char *)(y));
   endian_swap_long((unsigned char *)(z));
  }
 }

 ply_readOverhead(in, format, oh);

 return 1;
}

int ply_readFIndices(FILE *in, int format, int ph, int *nv, int *x, int *y, int *z)
{
 unsigned char nvs;
 int vc[3];

 ply_readOverhead(in, format, ph);

 if (format == PLY_FORMAT_ASCII) {fscanf(in,"%d %d %d %d", nv, x, y, z); return 1;}

 fread(&nvs, 1, 1, in);
 *nv = (int)nvs;
 fread(vc, 4, 3, in);
 *x = vc[0]; *y = vc[1]; *z = vc[2];

 if (format == PLY_FORMAT_BIN_B)
 {
  endian_swap_long((unsigned char *)(x));
  endian_swap_long((unsigned char *)(y));
  endian_swap_long((unsigned char *)(z));
 }

 return 1;
}

int ply_readAnotherFIndex(FILE *in, int format, int *x)
{
 if (format == PLY_FORMAT_ASCII) return (fscanf(in,"%d", x));

 if (fread(x, 4, 1, in) < 1) TMesh::error("Unexpected end of file!\n");

 if (format == PLY_FORMAT_BIN_B) endian_swap_long((unsigned char *)(x));

 return 1;
}

int Basic_TMesh::loadPLY(const char *fname)
{
 int format=0, voh, foh, vph, fph;
 int nv,nt,i,j,i1,i2,i3,i4;
 float x,y,z;
 bool triangulate = 0;
 FILE *in;
 char keyword[64], formats[24], version[10];
 Vertex *v;
 Node *n;

 if ((in = fopen(fname,"rb")) == NULL) TMesh::error("Can't open input ply file\n");

 if (strcmp(readLineFromFile(in),"ply")) TMesh::error("Input doesn't seem a valid ply file.\n");
 if (sscanf(readLineFromFile(in),"%7s %24s %10s",keyword,formats,version) < 3) TMesh::error("Unexpected token or end of file!\n");
 if (strcmp(keyword,"format")) TMesh::error("format definition expected!\n");
 if (!strcmp(formats,"ascii")) format = PLY_FORMAT_ASCII;
 else if (!strcmp(formats,"binary_little_endian")) format = PLY_FORMAT_BIN_L;
 else if (!strcmp(formats,"binary_big_endian")) format = PLY_FORMAT_BIN_B;
 else TMesh::error("Unrecognized format '%s'\n",formats);

 nv = ply_parseElements(in, "vertex");
 vph = ply_getOverhead(in, format, "vertex");
 ply_checkVertexProperties(in);
 voh = ply_getOverhead(in, format, "vertex");
 nt = ply_parseElements(in, "face");
 fph = ply_getOverhead(in, format, "face");
 ply_checkFaceProperties(in);
 foh = ply_getOverhead(in, format, "face");

 if (!sscanf(readLineFromFile(in),"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");
 while (strcmp(keyword, "end_header"))
  if (!sscanf(readLineFromFile(in),"%64s ",keyword)) TMesh::error("Unexpected token or end of file!\n");

 for (i=0; i<nv; i++)
 {
  ply_readVCoords(in, format, vph, voh, &x, &y, &z);
  V.appendTail(newVertex(x,y,z));
 }

 ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);
 i=0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

 i=0;
 TMesh::begin_progress();
 for (i=0; i<nt; i++)
 {
  if (ply_readFIndices(in, format, fph, &i4, &i1, &i2, &i3))
  {
   if ((i%1000) == 0) TMesh::report_progress("Loading ..%d%%",(i*100)/(nv*2));
   if (i1<0 || i2<0 || i3<0 || i4<3 || i1>(nv-1) || i2>(nv-1) || i3>(nv-1)) TMesh::error("\nloadPLY: Invalid index at face %d!\n",i);
   for (j=3; j<=i4; j++)
   {
    if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadPLY: Coincident indexes at triangle %d! Skipping.\n",i);
    else if (!CreateIndexedTriangle(var, i1, i2, i3)) TMesh::warning("\nloadPLY: This shouldn't happen!!! Skipping triangle.\n");
    i2 = i3;
    if (j<i4)
    {
     if (!ply_readAnotherFIndex(in, format, &i3)) TMesh::error("\nloadPLY: Couldn't read indexes for face # %d\n",i);
     else triangulate=1;
    }
    else ply_readOverhead(in, format, foh);
   }
  }
  else TMesh::error("\nloadPLY: Couldn't read indexes for face # %d\n",i);
 }
 TMesh::end_progress();

 closeLoadingSession(in, i, var, triangulate);
 TMesh::setFilename(fname);

 return 0;
}


////////////////////// Saves PLY format ///////////////////////////

int Basic_TMesh::savePLY(const char *fname, bool ascii)
{
 FILE *fp;
 int i, ii[3];
 float fc[3];
 char triname[256];
 unsigned char ii0 = 3;
 Node *n;
 coord *ocds;
 Vertex *v;

 strcpy(triname,fname);

 if ((fp = fopen(triname,"w")) == NULL)
 {
  TMesh::warning("Can't open '%s' for output !\n",triname);
  return 1;
 }

 fprintf(fp,"ply\n");
 if (ascii) fprintf(fp,"format ascii 1.0\n");
 else fprintf(fp,"format binary_little_endian 1.0\n");
 PRINT_PLY_COMMENT(fp);
 fprintf(fp,"element vertex %d\n",V.numels());
 fprintf(fp,"property float x\n");
 fprintf(fp,"property float y\n");
 fprintf(fp,"property float z\n");
 fprintf(fp,"element face %d\n",T.numels());
 fprintf(fp,"property list uchar int vertex_indices\n");
 fprintf(fp,"end_header\n");

 if (ascii) FOREACHVERTEX(v, n) fprintf(fp, "%f %f %f\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));
 else FOREACHVERTEX(v, n)
 {
  fc[0] = TMESH_TO_FLOAT(v->x); fc[1] = TMESH_TO_FLOAT(v->y); fc[2] = TMESH_TO_FLOAT(v->z);
  fwrite(fc, sizeof(float), 3, fp);
 }

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
 i=0; FOREACHVERTEX(v, n) v->x = i++;

 if (ascii) FOREACHNODE(T, n) fprintf(fp,"3 %d %d %d\n",TVI1(n),TVI2(n),TVI3(n));
 else FOREACHNODE(T, n)
 {
  ii[0]=TVI1(n); ii[1]=TVI2(n); ii[2]=TVI3(n);
  fwrite(&ii0, sizeof(unsigned char), 1, fp);
  fwrite(ii, sizeof(int), 3, fp);
 }

 fclose(fp);
 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}


////////////////////// Loads OBJ format ///////////////////////////

int Basic_TMesh::loadOBJ(const char *fname)
{
 FILE *fp;
 Node *n;
 char c, cmd[3] = "";
 float x,y,z;
 bool face_section = 0;
 int i=0,i1,i2,i3,nv=0,triangulate=0;
 Vertex *v;
 ExtVertex **var=NULL;

 if ((fp = fopen(fname,"r")) == NULL) return IO_CANTOPEN;

 TMesh::begin_progress();
 while (fscanf(fp, "%2s", cmd) && cmd[0] != '\0')
 {
  if (!strcmp(cmd,"v"))
  {
   if (face_section) TMesh::error("\nloadOBJ: Sorry. Couldn't manage disconnected vertex sections.\n");
   if (fscanf(fp, "%f %f %f", &x, &y, &z) == 3) V.appendTail(newVertex(x,y,z));
   else TMesh::error("\nloadOBJ: Couldn't read coordinates for vertex # %d\n",i);
  }
  else if (!strcmp(cmd,"f"))
  {
   if (!face_section)
   {
    nv = V.numels();
    var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);
    i=0; FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);
    face_section = 1;
    i=0;
   }

   if (fscanf(fp,"%d %d %d",&i1,&i2,&i3) == 3)
   {
    if ((i%1000) == 0) TMesh::report_progress("Loading ..%d%%",(i*100)/(nv*2));
    if (i1<0 || i2<0 || i3<0) TMesh::error("\nloadOBJ: Sorry. Negative vertex references not supported.\n");
    if (i1<1 || i2<1 || i3<1 || i1>nv || i2>nv || i3>nv) TMesh::error("\nloadOBJ: Invalid index at face %d!\n",i);
    do
    {
     if (i1 == i2 || i2 == i3 || i3 == i1) TMesh::warning("\nloadOBJ: Coincident indexes at triangle %d! Skipping.\n",i);
     else if (!CreateIndexedTriangle(var, i1-1, i2-1, i3-1)) TMesh::warning("\nloadOBJ: This shouldn't happen!!! Skipping triangle.\n");
     i2 = i3;
     while ((c=fgetc(fp)) != EOF && isspace(c) && c != '\n' && c != '\r');
     if (c==EOF) TMesh::error("\nloadOBJ: Unexpected end of file!\n");
     if (c != '\n' && c != '\r')
     {
      ungetc(c, fp);
      if (fscanf(fp,"%d",&i3) != 1) TMesh::error("\nloadOBJ: Couldn't read indexes for face # %d\n",i);
      else triangulate=1;
     }
    } while (c != '\n' && c != '\r');
   }
   else TMesh::error("\nloadOBJ: Couldn't read indexes for face # %d\n",i);
   i++;
  }
  else if (readLineFromFile(fp, 0) == NULL) break;
  cmd[0]='\0';
 }

 TMesh::end_progress();

 closeLoadingSession(fp, i, var, (triangulate!=0));
 TMesh::setFilename(fname);

 return 0;
}

////////////////////// Saves OBJ format ///////////////////////////

int Basic_TMesh::saveOBJ(const char *fname)
{
 FILE *fp;
 int i;
 char triname[256];
 Node *n;
 coord *ocds;
 Vertex *v;

 strcpy(triname,fname);

 if ((fp = fopen(triname,"w")) == NULL)
 {
  TMesh::warning("Can't open '%s' for output !\n",triname);
  return 1;
 }

 PRINT_HEADING_COMMENT(fp);

 FOREACHVERTEX(v, n) fprintf(fp, "v %f %f %f\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));

 ocds = new coord[V.numels()];
 i=0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
 i=0; FOREACHVERTEX(v, n) v->x = i++;

 FOREACHNODE(T, n) fprintf(fp,"f %d %d %d\n",TVI1(n)+1,TVI2(n)+1,TVI3(n)+1);

 fclose(fp);
 i=0; FOREACHVERTEX(v, n) v->x = ocds[i++];
 delete[] ocds;

 return 0;
}



////////////////////// Loads STL format ///////////////////////////

int Basic_TMesh::loadSTL(const char *fname)
{
 FILE *fp;
 int nt=0, i=0;
 char kw[64]="", kw2[64]="", *line, facet[50];
 float x,y,z;
 bool binary=0;
 Vertex *v, *v1=NULL, *v2=NULL, *v3=NULL;
 Edge *e1, *e2, *e3;
 Triangle *t;
 Point nor;

 if ((fp = fopen(fname,"r")) == NULL) return IO_CANTOPEN;

 fscanf(fp,"%5s",kw);
 if (strcmp(kw,"solid")) binary=1;

 if (!binary) // Might look ASCII just because it begins with 'solid'...
 {
	 rewind(fp);
	 if ((line = readLineFromFile(fp, 0)) == NULL) binary = 1;
	 else if ((line = readLineFromFile(fp, 0)) == NULL) binary = 1;
	 else {
		 sscanf(line, "%64s", kw);
		 if (strcmp(kw, "facet")) binary = 1;
	 }
 }
 if (binary) fp = freopen(fname, "rb", fp);

 TMesh::begin_progress();

 if (binary)
 {
  fseek(fp, 80, SEEK_SET);
  fread(&nt, 4, 1, fp);

  for (i=0; i<nt; i++)
  {
   if ((i)%10000 == 0) TMesh::report_progress(NULL);
   if (!fread(facet, 50, 1, fp)) TMesh::error("loadSTL: Unexpected end of file!\n");
   nor.setValue((*((float *)(facet+0))), (*((float *)(facet+4))), (*((float *)(facet+8))));
   v1 = newVertex((*((float *)(facet+12))), (*((float *)(facet+16))), (*((float *)(facet+20))));
   v2 = newVertex((*((float *)(facet+24))), (*((float *)(facet+28))), (*((float *)(facet+32))));
   v3 = newVertex((*((float *)(facet+36))), (*((float *)(facet+40))), (*((float *)(facet+44))));
   V.appendHead(v1); V.appendHead(v2); V.appendHead(v3);
   e1=CreateEdge(v1, v2); e2=CreateEdge(v2, v3); e3=CreateEdge(v3, v1);
   if (Triangle(e1,e2,e3).getNormal()*nor < 0) t=CreateTriangle(e1, e3, e2);
   else t=CreateTriangle(e1, e2, e3);
  }
 }
 else
 {
  while ((line = readLineFromFile(fp, 0))!=NULL)
  {
   if ((i++)%10000 == 0) TMesh::report_progress(NULL);
   sscanf(line,"%64s %f %f %f",kw,&x,&y,&z);
   if (!strcmp(kw,"facet"))
   {
    sscanf(line,"%64s %64s %f %f %f",kw,kw2,&x,&y,&z);
    nor.setValue(x,y,z);
   }
   else
   if (!strcmp(kw,"vertex"))
   {
    V.appendHead((v=newVertex(x,y,z)));
    if (v1==NULL) v1=v;
    else if (v2==NULL) v2=v;
    else if (v3==NULL)
    {
     v3=v;
     e1=CreateEdge(v1, v2);
     e2=CreateEdge(v2, v3);
     e3=CreateEdge(v3, v1);
     if (Triangle(e1,e2,e3).getNormal()*nor < 0) t=CreateTriangle(e1, e3, e2);
     else t=CreateTriangle(e1, e2, e3);
     v1=v2=v3=NULL;
    }
   }
  }
 }

 TMesh::end_progress();

 fclose(fp);

 TMesh::info("Loaded %d vertices and %d faces.\n",V.numels(),T.numels());

 if (T.numels()) rebuildConnectivity(); else return IO_UNKNOWN;
 TMesh::setFilename(fname);

 return 0;
}

////////////////////// Saves STL format ///////////////////////////

int Basic_TMesh::saveSTL(const char *fname)
{
	FILE *fp;
	char triname[256];
	Node *n;
	Triangle *t;
	Point nor;

	strcpy(triname, fname);

	if ((fp = fopen(triname, "w")) == NULL)
	{
		TMesh::warning("Can't open '%s' for output !\n", triname);
		return 1;
	}

	fprintf(fp, "solid T_MESH\n");

	FOREACHTRIANGLE(t, n)
	{
		nor = t->getNormal();
		fprintf(fp, " facet normal %f %f %f\n", TMESH_TO_FLOAT(nor.x), TMESH_TO_FLOAT(nor.y), TMESH_TO_FLOAT(nor.z));
		fprintf(fp, "  outer loop\n");
		fprintf(fp, "   vertex %f %f %f\n", TMESH_TO_FLOAT(t->v1()->x), TMESH_TO_FLOAT(t->v1()->y), TMESH_TO_FLOAT(t->v1()->z));
		fprintf(fp, "   vertex %f %f %f\n", TMESH_TO_FLOAT(t->v2()->x), TMESH_TO_FLOAT(t->v2()->y), TMESH_TO_FLOAT(t->v2()->z));
		fprintf(fp, "   vertex %f %f %f\n", TMESH_TO_FLOAT(t->v3()->x), TMESH_TO_FLOAT(t->v3()->y), TMESH_TO_FLOAT(t->v3()->z));
		fprintf(fp, "  endloop\n");
		fprintf(fp, " endfacet\n");
	}
	fprintf(fp, "endsolid T_MESH\n");

	fclose(fp);

	return 0;
}

////////////////////// Saves Extact File format (EFF) ///////////////////////////

int Basic_TMesh::saveEFF(const char *fname)
{
	std::ofstream os;
	os.open(fname);
	if (!os.is_open())
	{
		TMesh::warning("Can't open '%s' for output !\n", fname);
		return 1;
	}

	os << "EFF\n";
	os << V.numels() << " " << T.numels() << "\n";

	Node *n;
	Vertex *v;
	FOREACHVERTEX(v, n) os << v->x << " " << v->y << " " << v->z << "\n";

	coord *ocds = new coord[V.numels()];
	int i;
	i = 0; FOREACHVERTEX(v, n) ocds[i++] = v->x;
	i = 0; FOREACHVERTEX(v, n) v->x = i++;

	FOREACHNODE(T, n) os << TVI1(n) << " " << TVI2(n) << " " << TVI3(n) << "\n";

	os.close();
	i = 0; FOREACHVERTEX(v, n) v->x = ocds[i++];
	delete[] ocds;

	return 0;
}

} //namespace T_MESH
