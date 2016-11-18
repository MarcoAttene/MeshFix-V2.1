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

#ifndef _TIN_H
#define _TIN_H

#include "tmesh.h"

namespace T_MESH
{

//! Basic_TMesh

//! This class represents a manifold and oriented triangle mesh.
//! Vertices, Edges and Triangles are stored in the Lists V, E and T
//! respectively. Methods boundaries(), handles() and shells() may
//! be used to retrieve the respective topological entities.
//! Navigation of the mesh is based on the topological relationships
//! stored in each Vertex, Edge and Triangle.
//! Some methods would require a global update only to maintain
//! consistent values of the protected fields n_boundaries, n_handles
//! and n_shells. On the other hand, the same methods would work in
//! constant time if these values would not need to be updated.
//! To keep complexity as low as possible,
//! we make use of the 'dirty bits' d_boundaries, d_handles and 
//! d_shells to mark that the respective entities must be updated.
//! The access functions boundaries(), handles() and shells()
//! check the status of the respective dirty bit and do a global
//! update (i.e., eulerUpdate()) only if necessary.
//! The complexity of the methods is provided as a function of a
//! generic 'N', which is O(V.numels()) = O(E.numels()) = O(T.numels()).

class Basic_TMesh
{
	protected:

		int n_boundaries;		//!< Number of boundary loops
		int n_handles;			//!< Number of handles
		int n_shells;			//!< Number of connected components

		bool d_boundaries;		//!< Dirty bit for n_boundaries
		bool d_handles;		//!< Dirty bit for n_handles
		bool d_shells;			//!< Dirty bit for n_shells

	public:

		List V;			//!< Vertex set
		List E;			//!< Edge set
		List T;			//!< Triangle set

		void *info;		//! Generic information attached to this mesh

		/////////////////////////////////////////////////////////////////////////////
		//
		// Constructors/Destructor (Implemented in "MESH_STRUCTURE/io.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Empty triangulation. Should be used only prior to a call to load().
		Basic_TMesh();

		//! Pre-defined triangulation. Currently, only "triangle" and "tetrahedron"
		//! are recognized.
		Basic_TMesh(const char *);
		void init(const char *);

		//! Clones an existing Trianglation.
		Basic_TMesh(const Basic_TMesh *, const bool clone_info = false);
		void init(const Basic_TMesh *, const bool clone_info = false);

		//! Clones an existing connected component.

		//! Creates a new Basic_TMesh out of a connected component of an existing
		//! Basic_TMesh. 't' is a triangle of the connected component that must
		//! be copied. If 'keep_ref' is TRUE, each element of the existing mesh
		//! keeps a pointer to the corresponding new element in the 'info' field.
		Basic_TMesh(const Triangle *t, const bool keep_ref = false);
		void init(const Triangle *t, const bool keep_ref = false);

		//! Destructor. Frees the memory allocated for all the mesh elements.
		//! Warning! This method uses the freeNodes() method of the class List,
		//! which is not guaranteed to work correctly on systems other than
		//! Linux (see the documentation of List for details).
		//! Assuming that the method works correctly, however, the calling
		//! function is responsible of freeing the memory that was possibly
		//! allocated for objects pointed to by the 'info' field of mesh
		//! elements. Clearly, this must be done before calling the destructor.
		~Basic_TMesh();

		//! Returns true only if object is a basic Basic_TMesh. All the reimplementations must return false.
//		TMESH_VIRTUAL bool isBaseType() const { return true; }

		//! Get the number of boundary loops of the triangle mesh. O(1) or O(N).
		int boundaries() { if (d_boundaries) eulerUpdate(); return n_boundaries; }

		//! Get the number of handles of the triangle mesh. O(1) or O(N).
		int handles() { if (d_handles) eulerUpdate(); return n_handles; }

		//! Get the number of connected components of the triangle mesh. O(1) or O(N).
		int shells() { if (d_shells) eulerUpdate(); return n_shells; }


		/////////////////////////////////////////////////////////////////////////////
		//
		// Input/Output methods (Implemented in "MESH_STRUCTURE/io.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Initialize the triangle mesh from the file 'filename'.

		//! The file format is automatically deduced from the magic number
		//! or the filename extension. If 'doupdate' is FALSE, the
		//! global update for the topological entities is prevented.
		//! Currently, the following file formats are supported:
		//! Open Inventor (IV), VRML 1.0 and 2.0 (WRL), Object File Format (OFF),
		//! IMATI Ver-Tri (VER, TRI), PLY, OBJ, STL.
		//! A non-zero value is returned in case of error. Specifically,
		//! IO_CANTOPEN means that the file couldn't be opened for reading.
		//! IO_FORMAT means that the file format was not recognized by the loader.
		//! IO_UNKNOWN represents all the other errors.
		//! The calling function is responsible of verifying that the mesh is
		//! empty before calling this method.

		int load(const char *filename, const bool update = 1);
		int loadIV(const char *);		//!< Loads IV
		int loadVRML1(const char *);		//!< Loads VRML 1.0
		int loadOFF(const char *);		//!< Loads OFF
		int loadEFF(const char *);		//!< Loads EFF
		int loadPLY(const char *);		//!< Loads PLY
		int loadVerTri(const char *);		//!< Loads VER-TRI
		int loadVRML2(const char *);		//!< Loads VRML 2.0
		int loadOBJ(const char *);		//!< Loads OBJ
		int loadSTL(const char *);		//!< Loads STL

		int cutAndStitch();	//!< Convert to manifold
		Triangle * CreateIndexedTriangle(ExtVertex **, int, int, int);							
		TMESH_VIRTUAL Triangle * CreateTriangleFromVertices(ExtVertex *, ExtVertex *, ExtVertex *);	

		//! This function approximates the vertex coordinates with the values
		//! that can be represented in an ASCII file.
		void coordBackApproximation();

	protected:
		void closeLoadingSession(FILE *, int, ExtVertex **, bool);
		bool pinch(Edge *e, bool wcv);
		//! If the 'e' is internal, creates a copy of the edge and associates e->t2
		//! to this copy. After this operation, 'e' and its copy are boundary edges
		//! with the same vertices. If 'e' is already on boundary, nothing is done
		//! and NULL is returned. Otherwise the newly created copy is returned.
		TMESH_VIRTUAL Edge *duplicateEdge(Edge *e);

	public:
		TMESH_VIRTUAL Vertex *	newVertex();
		TMESH_VIRTUAL Vertex *	newVertex(const coord &, const coord &, const coord &);
		TMESH_VIRTUAL Vertex *	newVertex(Point *);
		TMESH_VIRTUAL Vertex *	newVertex(Point &);
		TMESH_VIRTUAL Vertex *	newVertex(Vertex *);
		TMESH_VIRTUAL Edge *		newEdge(Vertex *, Vertex *);
		TMESH_VIRTUAL Edge *		newEdge(Edge *);
		TMESH_VIRTUAL Triangle *	newTriangle();
		TMESH_VIRTUAL Triangle *	newTriangle(Edge *, Edge *, Edge *);
		TMESH_VIRTUAL Basic_TMesh *	newObject() const { return new Basic_TMesh(); }
		TMESH_VIRTUAL Basic_TMesh *	newObject(const Basic_TMesh *tm, const bool ci = false) const { return new Basic_TMesh(tm, ci); }
		TMESH_VIRTUAL Basic_TMesh *   newObject(const char *s) const { return new Basic_TMesh(s); }

		//! Save the triangle mesh to file 'filename'.

		//! The file format is deduced from the filename extension
		//! (wrl = vrml 1.0), (iv = OpenInventor), (off = Object
		//! file format), (ply = PLY format), (tri = IMATI Ver-Tri).
		//! If 'back_approx' is set, vertex coordinates are approximated
		//! to reflect the limited precision of floating point
		//! representation in ASCII files. This should be used when
		//! coherence is necessary between in-memory and saved data.
		//! A non-zero return value is returned if errors occur.

		int save(const char *filename, bool back_approx = 0);

		int saveIV(const char *);		//!< Saves IV
		int saveOFF(const char *);		//!< Saves OFF 1.0
		int saveEFF(const char *);		//!< Saves EFF
		int saveOBJ(const char *);		//!< Saves OBJ
		int saveSTL(const char *);		//!< Saves STL
		int savePLY(const char *, bool ascii = 1); //!< Saves PLY 1.0 (ascii or binary)
		int saveVerTri(const char *);		//!< Saves Ver-Tri

		//! Saves the triangle mesh to a VRML 1.0 file.
		//! The value of 'mode' specifies whether to use additional
		//! information attached to mesh elements in order to assign
		//! them a proper color.
		//! IO_CSAVE_OVERALL assigns a unique color for the entire mesh (default).
		//! IO_CSAVE_PERFACE assigns a color to each triangle depending on the value
		//! of its 'info' field.
		//! IO_CSAVE_PERVERTEX	assigns a color to each vertex depending on the value
		//! of its 'info' field.
		//! IO_CSAVE_PERFACE_INDEXED assigns one of five base colors to each triangle
		//! depending on the value of its 'mask' field.
		//! IO_CSAVE_PERVERTEX_INDEXED assigns one of five base colors to each vertex
		//! depending on the value of its 'mask' field.
		int saveVRML1(const char *, const int mode = 0);


		//! Append another triangle mesh to the existing one.

		//! This method works exactly as the 'load()' method, except for the fact
		//! that it does not assume that the mesh is empty.
		int append(const char *filename, const bool doupdate = 1);


		// Move all the elements of 't' to this mesh and delete 't' itself.
		TMESH_VIRTUAL void moveMeshElements(Basic_TMesh *t, bool delInput = true);	


		/////////////////////////////////////////////////////////////////////////////
		//
		// Primitive Construction (Implemented in "MESH_STRUCTURE/tin.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Creates an Edge connecting two existing mesh vertices.

		//! Returns the newly created edge. If an edge connecting the two vertices
		//! already exists in the mesh, then no new edge is created and the old one
		//! is returned.
		Edge     *CreateEdge(Vertex *v1, Vertex *v2);


		//! Creates an Edge connecting two existing mesh Extended vertices.

		//! Returns the newly created edge. If an edge connecting the two vertices
		//! already exists in the mesh, then no new edge is created and the old one
		//! is returned.
		//! If 'check' is FALSE, the check for previously existing edges is skipped.
		Edge     *CreateEdge(ExtVertex *v1, ExtVertex *v2, const bool check = 1);


		//! Creates a properly oriented Triangle bounded by three existing mesh edges.

		//! Returns the newly created Triangle. If e1, e2 and e3
		//! are not suitable for creating a properly oriented and
		//! manifold triangle, the creation fails and NULL is returned.
		TMESH_VIRTUAL Triangle * CreateTriangle(Edge *e1, Edge *e2, Edge *e3);		


		//! Creates an arbitrarily oriented Triangle bounded by three existing mesh edges.

		//! Returns the newly created Triangle. If either e1, e2 or e3
		//! has already two incident triangles, the creation fails and NULL is returned.
		//! This method assumes that e1, e2 and e3 are incident to exactly three vertices.
		TMESH_VIRTUAL Triangle * CreateUnorientedTriangle(Edge *, Edge *, Edge *);	


		//! Creates a newEdge 'e' and an oriented Triangle bounded by 'e', 'e1' and 'e2'.

		//! The newly created triangle is returned, unless 'e1' and 'e2' do not share a
		//! vertex or they are not boundary edges. In this cases, NULL is returned.
		TMESH_VIRTUAL Triangle *EulerEdgeTriangle(Edge *e1, Edge *e2);			


		//! Splits and edge at a given point and returns the newly created Vertex.
		//! If the boolean parameter is set to true, the 'mask' fields of edges and
		//! triangles are propagated to the new elements.
		TMESH_VIRTUAL Vertex *splitEdge(Edge *, Point *, bool = 0);


		//! Splits a triangle at a given point and returns the newly created Vertex.
		//! If the boolean parameter is set to true, the 'mask' field of the
		//! triangle is propagated to the new triangle.
		TMESH_VIRTUAL Vertex *splitTriangle(Triangle *, Point *, bool = 0);

		//! Creates two new triangles connecting the boundary edges e1 and e2
		//! and returns their common edge.
		//! If e1 and e2 share a vertex, then only one triangle is created and
		//! e1 is returned.
		//! Returns NULL if either e1 or e2 are not boundary edges.
		TMESH_VIRTUAL Edge *bridgeBoundaries(Edge *e1, Edge *e2);					

		/////////////////////////////////////////////////////////////////////////////
		//
		// Primitive Destruction (Implemented in "MESH_STRUCTURE/tin.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////


		//! Unlinks a triangle from the mesh. O(1).

		//! Resulting isolated vertices and edges are unlinked too.
		//! If necessary, this method duplicates non-manifold vertices that can
		//! occur due to the removal of the triangle.
		//! The unlinked triangle, along with the other possible unlinked elements,
		//! must be removed from the List T through removeUnlinkedElements().
		void unlinkTriangle(Triangle *);


		//! Unlinks a triangle from the mesh. O(1).

		//! No check is performed on the resulting topology, which may be inconsistent.
		//! The unlinked triangle, along with the other possible unlinked elements,
		//! must be removed from the List T through removeUnlinkedElements().
		void unlinkTriangleNoManifold(Triangle *);


		//! Removes a triangle from the mesh. O(N).

		//! This is equivalent to an unlinkTriangle(t) followed by a
		//! removeUnlinkedElements().
		void removeTriangle(Triangle *t) { unlinkTriangle(t); removeUnlinkedElements(); }

		//! Removes all the unlinked triangles from List T. Returns the number of removed triangles. O(N).
		int removeTriangles();

		//! Removes all the unlinked edges from List E. Returns the number of removed edges. O(N).
		int removeEdges();

		//! Removes all the unlinked vertices from List V. Returns the number of removed vertices. O(N).
		int removeVertices();

		//! Removes all the unlinked elements from the lists. Returns the number of removed elements. O(N).
		int removeUnlinkedElements() { return removeTriangles() + removeEdges() + removeVertices(); }

		//! Removes all the vertices that can be deleted without changing the geometric realization. O(N).
		int removeRedundantVertices();

		/////////////////////////////////////////////////////////////////////////////
		//
		// Methods acting on selections (Implemented in "MESH_STRUCTURE/tin.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Deselects all the triangles. O(N).
		void deselectTriangles();

		//! Removes all the selected triangles. O(N).
		void removeSelectedTriangles();

		//! Selects all the triangles having at least one boundary vertex. O(N).
		//! Returns the number of selected triangles.
		int selectBoundaryTriangles();

		//! Enlarges the current selection of one triangle in width. O(N).

		//! Each triangle sharing at least one vertex with a currently selected
		//! triangle becomes selected.
		//! Returns the number of newly selected triangles.
		int growSelection();

		//! Shrinks the current selection of one triangle in width. O(N).

		//! Each triangle sharing at least one vertex with a currently unselected
		//! triangle becomes unselected.
		void shrinkSelection();

		//! Inverts the selection status of all the triangles. O(N).

		//! If 't0' is not NULL, then only the connected component containing 't0'
		//! is inverted.
		void invertSelection(Triangle *t0 = NULL);

		//! If 't0' is selected, deselect everything but the selected triangles connected to 't0'
		void reselectSelection(Triangle *t0);

		//! Creates a new Basic_TMesh out of an existing selection containing 't0'. O(output).

		//! If necessary, non-manifold vertices are properly duplicated.
		//! If 'keep_ref' is set to TRUE, then elements of the original mesh point
		//! (through their info field) to corresponding elements of the newly created copy.

		TMESH_VIRTUAL Basic_TMesh *createSubMeshFromSelection(Triangle *t0 = NULL, bool keep_ref = 0);	

		//! Creates a new Basic_TMesh out of an existing triangle 't0'. O(output).

		TMESH_VIRTUAL Basic_TMesh *createSubMeshFromTriangle(Triangle *t0);	

		//! Marks all the triangles within distance L from 'p' as selected. O(output).

		//! A triangle is considered to be within distance L from 'p' only if all
		//! its three vertices are so.
		//! Point 'p' is assumed to belong to triangle 't0', which is required to
		//! limit the complexity to the size of the selected region.
		//! Returns the number of selected triangles.
		int  selectSphericalRegion(Triangle *t0, const double L, const Point *p);


		//! Marks all the triangles within distance L from 'p' as deselected. O(output).

		//! A triangle is considered to be within distance L from 'p' only if all
		//! its three vertices are so.
		//! Point 'p' is assumed to belong to triangle 't0', which is required to
		//! limit the complexity to the size of the selected region.
		//! Returns the number of deselected triangles.
		int  deselectSphericalRegion(Triangle *t0, const double L, const Point *p);


		//! Deselects all the triangles farther than L from 'p'. O(N).

		//! A triangle is considered to be farther than L from 'p' if at least
		//! one of its three vertices is so.
		//! Point 'p' is assumed to belong to triangle 't0'. Passing 't0' avoids
		//! the non robust and expensive computation of point-in-triangle.
		void reselectSphericalRegion(Triangle *t0, const double L, const Point *p);

		//! Re-triangulates the currently selected region using a Delaunay-like approach. O(SlogS).

		//! A common plane is computed as the average of the planes of the triangles selected;
		//! then, the vertices of the region are projected on the plane and edges are iteratively
		//! swapped up to convergence (which is guaranteed on planar and simple domains).
		//! Finally, the vertices are moved back to their original positions. This operation is
		//! particularly useful to improve the quality of nearly flat regions. The selection must
		//! be simple and its Gauss map must be entirely contained in a semi-sphere.
		//! Returns TRUE on success, FALSE otherwise.
		bool retriangulateSelectedRegion();


		//! TRUE iff the set of selected triangles in 'l' is simply connected. O(l->numels()).
		bool isSelectionSimple(List *l);

		//! Unmarks all the elements but leaves the selection status of triangles as is. O(N).
		void unmarkEverythingButSelections();


		//! Selects all the triangles of the connected component containing t0. O(N).

		//! If 'stop_on_sharp', expansion from 't0' brakes at tagged sharp edges.
		//! Returns the number of selected triangles.
		int selectConnectedComponent(Triangle *t0, bool stop_on_sharp = 0);


		//! Deselects all the triangles of the connected component containing t0. O(N).

		//! If 'stop_on_sharp', expansion from 't0' brakes at tagged sharp edges.
		//! Returns the number of deselected triangles.
		int deselectConnectedComponent(Triangle *t0, bool stop_on_sharp = 0);


		//! Append to the current mesh a copy of all the elements of 't'.
		//! The newly created elements form a new selection.
		TMESH_VIRTUAL void append(Basic_TMesh *t);	

		//! This method removes one connected component from the mesh and creates
		//! a separate new mesh out of it. The components to be removed is the one
		//! containing the first triangle in the list T.
		//! Possible selection flags are deleted by this method.
		TMESH_VIRTUAL Basic_TMesh *split();			


		/////////////////////////////////////////////////////////////////////////////
		//
		// Region manipulation (Implemented in "MESH_STRUCTURE/tin.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Make a list of triangles within distance L from 'p'. O(output).

		//! Starting from 't0', which is assumed to contain 'p', add a triangle at a
		//! time to the list as long as all the vertices stay within distance L from 'p'.
		List     *getRegion(Triangle *t0, const double L, const Point *p);

		//! Removes triangles within distance L from 'p'. O(N).

		//! Starting from 't0', which is assumed to contain 'p', remove a triangle at a
		//! time as long as all its vertices stay within distance L from 'p'.
		void      removeRegion(Triangle *t0, const double L, const Point *p);

		//! Get the vertex next to 'v' on the boundary of the region. O(1).
		Vertex   *nextVertexOnRegionBoundary(Vertex *v) const;

		//! Retrieve internal vertices of a region. O(l->numels()).

		//! This method returns a list containing an edge of the region's boundary
		//! as its first element, and all the internal vertices as the remaining elements.
		List     *getRegionInternalVertices(List *l);

		//! Transform the vertices of the shell containing 't0' using the matrix m. O(S).
		void      transformShell(Triangle *t0, const Matrix4x4& m);

		//! Translate the mesh by a vector
		void translate(const Point& t_vec);

		//! Return the center of mass of the mesh
		Point getCenter() const;

		//! Remove all the triangles belonging to the shell containing 't0'. O(N).
		void 	   removeShell(Triangle *t0);


		/////////////////////////////////////////////////////////////////////////////
		//
		// Global Operations (Implemented in "MESH_STRUCTURE/tin.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Tag sharp edges based on threshold angle. O(N).

		//! Tag as sharp all the edges in which the normals of the two incident
		//! triangles form an angle greater than 't'.
		void   sharpEdgeTagging(const double t);

		//! Unmark all the elements. O(N).
		void   unmarkEverything();

		//! Bounding box longest edge. 'b' and 't' are set as the longest diagonal end-points. O(N).
		coord getBoundingBox(Point& b, Point& t) const;

		//! Bounding box longest diagonal. O(N).
		double bboxLongestDiagonal() { Point a, b; getBoundingBox(a, b); return a.distance(b); }

		//! Approximate bounding ball radius. O(N).
		double getBoundingBallRadius() const;

		//! Total area of the mesh. O(N).
		double area() const;

		//! Total volume of the mesh assuming that boundaries() = 0. O(N).
		double volume() const;

		//! Scale the mesh to make it fit within a cube [0,0,0]-[s,s,s]. O(N).
		void   normalize(const coord s = 1.0);

		//! Scale the mesh to make it fit within a cube [0,0,0]-[s,s,s] and snap coordinates on grid points O(N).
		void quantize(const int s = 65536);

		//! Transform the mesh geometry using the transformation matrix m. O(N).
		void   transform(const Matrix4x4& m);

		//! Randomly move vertices along their normals. O(N).
		//! Displacement is bounded by 'p'% of the bounding ball radius.
		void   addNormalNoise(const double p);

		//! Iteratively swaps edges to minimize the Delaunay minimum angle. O(N).

		//! Edges tagged as sharp are constrained not to swap.
		//! On generically curved manifolds this process is not guaranteed to converge.
		//! This method returns TRUE if convergence is reached, FALSE otherwise.
		bool   iterativeEdgeSwaps();

		//! True if the mesh properly contains 'p' (exact when using rationals).
		//! The mesh is assumed to be a well-defined polyhedron (i.e. no boundary,
		//! no degenerate triangles, no self-intersections) and to have a correct orientation:
		//! result is undetermined otherwise.
		bool isInnerPoint(Point& p) const;

		//! Performs one step of Loop subdivision. If 'midpoint' is set, only the connectivity
		//! is subdivided, whereas the surface shape is kept unchanged.
		void loopSubdivision(bool midpoint =false);

		/////////////////////////////////////////////////////////////////////////////
		//
		// Surface topology manipulation (Implemented in "MESH_STRUCTURE/tin.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Invert all the triangle and edge orientations. O(N).
		void      flipNormals();

		//! Invert the orientation of triangles and edges belonging to the shell containing 't0'. O(S).
		void      flipNormals(Triangle *t0);

		//! Return the triangle with the maximum 'z' coordinate in the shell containing 't0'. O(N).

		//! Useful for orienting meshes bounding solids.
		Triangle *topTriangle(Triangle *t0);

		//! Updates the values of n_boundaries, n_handles and n_shells. O(N).

		//! The relative dirty bits are set to zero.
		void      eulerUpdate();

		//! Duplicates edges and vertices to make the mesh homeomorphic to a disk. O(N).
		void      openToDisk();


		/////////////////////////////////////////////////////////////////////////////
		//
		// Topological fixes (Implemented in "MESH_STRUCTURE/checkAndRepair.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Removes all connected components but the one having most triangles.
		//! Returns the number of components removed.
		int       removeSmallestComponents();

		//! Checks that triangles are consistently oriented and, if they are not,
		//! invert some of them to achieve an overall consistency. If the mesh is
		//! not orientable cut it. Returns: 0 if mesh was already oriented; 1 if
		//! the mesh could be oriented without cuts; >1 if cuts were necessary.
		int       forceNormalConsistence();

		//! Same as above, but acts on a single connected component and uses one
		//! specific triangle from which the orientation is propagated.
		int       forceNormalConsistence(Triangle *);

		//! Detect singular vertices and duplicte them. Return number of singular
		//! vertices being duplicated.
		int       duplicateNonManifoldVertices();

		//! Remove redundant triangles (i.e. having the same vertices as others)
		//! and return their number.
		int       removeDuplicatedTriangles();

		//! Check the mesh connectivity. If everything is fine NULL is returned,
		//! otherwise an error string is returned.
		const char *checkConnectivity();

		//! If called in rebuildConnectivity(bool) fix the connectivity between 
		//! geometric elements
		bool fixConnectivity(); //!< AMF_ADD 1.1>

		//! Considers triangles as purely geometric entities and recomputes their
		//! connectivity based on vertex correspondence.
		//! Returns false if mesh was not an oriented manifold.
		TMESH_VIRTUAL bool rebuildConnectivity(bool = true); //!< AMF_CHANGE 1.1>

		//! Looks for topologically different edges having the same geometry
		//! (i.e. coincident vertices) and unify them. Return their number.
		int     mergeCoincidentEdges();


		/////////////////////////////////////////////////////////////////////////////
		//
		// Geometrical fixes (Implemented in "MESH_STRUCTURE/checkAndRepair.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Looks for zero-area triangles and resolve them by splits and collapses.
		//! If, for any reason, some degeneracies cannot be removed
		//! a warning is issued and the degenerate triangles that could not
		//! be resolved are selected.
		//! The absolute value of the integer returned is the number of 
		//! collapses performed; the return value is negative if some 
		//! degenerate triangles could not be resolved.
		int removeDegenerateTriangles();

		//! Calls 'removeDegenerateTriangles()' and, if some degeneracies remain,
		//! removes them and fills the resulting holes. Then tries again and, if
		//! some degeneracies still remain, removes them and their neighbors and
		//! fills the resulting holes, and so on, until the neighborhood growth
		//! reaches max-iters. If even in this case some degeneracies remain,
		//! returns false, otherwise returns true.
		bool strongDegeneracyRemoval(int max_iters);

		//! Removes all the intersecting triangles and patches the resulting holes.
		//! If the patches still produce intersections, iterates again on a larger
		//! neighborhood. Tries up to max_iters times before giving up. Returns
		//! true only if all the intersections could be removed.
		bool strongIntersectionRemoval(int max_iters);

		//! Iteratively call strongDegeneracyRemoval and strongIntersectionRemoval
		//! to produce an eventually clean mesh without degeneracies and intersections.
		//! The two aforementioned methods are called up to max_iter times and
		//! each of them is called using 'inner_loops' as a parameter.
		//! Returns true only if the mesh could be completely cleaned.
		bool meshclean(int max_iters = 10, int inner_loops = 3);

		//! Removes overlapping triangles and return their number.
		int removeOverlappingTriangles();

		//! Checks the mesh for degeneracies, concident vertices and overlaps.
		//! If such a flaw is found returns its closest vertex.
		Vertex *checkGeometry();

		//! Selects all the triangles that unproperly intersect other parts of
		//! the mesh and return their number. The parameter 'tris_per_cell'
		//! determines the depth of the recursive space subdivision used to keep
		//! the complexity under a resonable threchold. The default value is safe
		//! in most cases.
		//! if 'justproper' is true, coincident edges and vertices are not regarded
		//! as intersections even if they are not common subsimplexes.
		int selectIntersectingTriangles(UINT16 tris_per_cell = 50, bool justproper = false);


		//! This is as coordBackApproximation() but it also checks for
		//! intersections and, if any, it tries different approximations.
		//! Returns true if no intersections remain.
		bool safeCoordBackApproximation();

		//! Removes all the connected components whose area is less than 'epsilon'.
		//! Returns the number of components removed.
		int       removeSmallestComponents(double epsilon);


		/////////////////////////////////////////////////////////////////////////////
		//
		// Hole triangulation (Implemented in "MESH_STRUCTURE/holeFilling.cpp")
		//
		/////////////////////////////////////////////////////////////////////////////

		//! Computes the barycenter of the boundary loop and connects it with
		//! all the edges of the loop. Returns the number of triangles
		//! created.
		TMESH_VIRTUAL int StarTriangulateHole(Edge *);	

		//! Creates a triangulation whose projection on the plane with normal 'nor'
		//! is Delaunay. Returns the number of triangles created.
		int TriangulateHole(Edge *, Point *nor);

		//! Creates a triangulation of the hole based on heuristics.
		//! Returns the number of triangles created.
		int TriangulateHole(Edge *);

		//! Creates a triangulation of the hole based on the assumption that it is flat.
		//! Returns the number of triangles created.
		//		int TriangulateFlatHole(Edge *);

		//! Creates a triangulation and inserts the additional points in the List.
		//! Returns the number of triangles created.
		TMESH_VIRTUAL int TriangulateHole(Edge *, List *);			

		//! Hole filling algorithm. Performs a triangulation based on heuristics and,
		//! if 'refine' is set to true, adds inner vertices to reproduce the sampling
		//! density of the surroundings.
		void FillHole(Edge *, bool refine = true);

		//! Fills all the holes having at least 'nbe' boundary edges. If 'refine'
		//! is true, adds inner vertices to reproduce the sampling density
		//! of the surroundings. Returns number of holes patched.
		//! If 'nbe' is 0 (default), all the holes are patched.
		int fillSmallBoundaries(int nbe = 0, bool refine = true);

		//! Takes a selected region and inserts inner vertices to reproduce the
		//! sampling density of the surroundings. If 't0' is not NULL, only the
		//! selected region containing 't0' is refined. Returns the number of
		//! vertices inserted.
		TMESH_VIRTUAL int refineSelectedHolePatches(Triangle *t0 =NULL);		

		//! Retriangulates the vertex neghborhood based on heuristics.
		int retriangulateVT(Vertex *);

		//! Joins the two boundary vertices gv and gw through an edge. A pair of triangles is
		//! added to properly change the topology of the mesh.
		Edge	*joinBoundaryLoops(Vertex *, Vertex *, bool = 0, bool = 1);

		// Debug and work-in-progress

		void printReport();

	protected:
		Vertex *watsonInsert(Point *, List *, int);
	};

#define FOREACHTRIANGLE(Tt, n) for (n = T.head(), Tt = (n)?((Triangle *)n->data):NULL; n != NULL; n=n->next(), Tt = (n)?((Triangle *)n->data):NULL)
#define FOREACHEDGE(Tt, n) for (n = E.head(), Tt = (n)?((Edge *)n->data):NULL; n != NULL; n=n->next(), Tt = (n)?((Edge *)n->data):NULL)
#define FOREACHVERTEX(Tt, n) for (n = V.head(), Tt = (n)?((Vertex *)n->data):NULL; n != NULL; n=n->next(), Tt = (n)?((Vertex *)n->data):NULL)

#define MARK_VISIT(a)   ((a)->mask |= ((unsigned char)1))
#define IS_VISITED(a)   ((a)->mask &  ((unsigned char)1))
#define UNMARK_VISIT(a) ((a)->mask &= (~((unsigned char)1)))

#define MARK_VISIT2(a)   ((a)->mask |= ((unsigned char)2))
#define IS_VISITED2(a)   ((a)->mask &  ((unsigned char)2))
#define UNMARK_VISIT2(a) ((a)->mask &= (~((unsigned char)2)))

#define MARK_BIT(a,b)   ((a)->mask |= ((unsigned char)(1<<b)))
#define IS_BIT(a,b)     ((a)->mask &  ((unsigned char)(1<<b)))
#define UNMARK_BIT(a,b) ((a)->mask &= (~((unsigned char)(1<<b))))

#define TAG_SHARPEDGE(a)   (MARK_BIT((a),7))
#define IS_SHARPEDGE(a)    (IS_BIT((a),7))
#define UNTAG_SHARPEDGE(a) (UNMARK_BIT((a),7))


	// Errors from loading

#define IO_CANTOPEN	10
#define IO_FORMAT	20
#define IO_UNKNOWN	30

#define IO_CSAVE_OVERALL		0
#define IO_CSAVE_PERFACE		1
#define IO_CSAVE_PERVERTEX		2
#define IO_CSAVE_PERFACE_INDEXED	3
#define IO_CSAVE_PERVERTEX_INDEXED	4

} //namespace T_MESH

#endif //_TIN_H

