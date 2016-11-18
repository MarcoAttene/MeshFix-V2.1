----------------------------
MeshFIX - Version 2.1                                                          
----------------------------

by Marco Attene
                                                                    
Consiglio Nazionale delle Ricerche                                        
Istituto di Matematica Applicata e Tecnologie Informatiche                
Sezione di Genova                                                         
IMATI-GE / CNR                                                            
                                                                         
This software takes as input a polygon mesh and produces a copy of the input where all the occurrences of a specific set of "defects" are corrected. MeshFix has been designed to correct typical flaws present in RAW DIGITIZED mesh models, thus it might fail or produce coarse results if run on other sorts of input meshes (e.g. tessellated CAD models). When the software fails, it appends a textual description to the file "meshfix.log".

The input is assumed to represent a single CLOSED SOLID OBJECT, thus the output will be a SINGLE WATERTIGHT TRIANGLE MESH bounding a polyhedron. All the singularities, self-intersections and degenerate elements are removed from the input, while regions of the surface without defects are left unmodified. Accepted input file formats are:
STL (http://www.sdsc.edu/tmf/Stl-specs/stl.html)
OFF (http://shape.cs.princeton.edu/benchmark/documentation/off_format.html)
PLY (http://www.cs.unc.edu/~geom/Powerplant/Ply.doc)
and partially:
IV 2.1, VRML 1.0, VRML 2.0, OBJ.


-----------------------------
ALGORITHM AND CITATION POLICY
-----------------------------
To better understand how the algorithm works, please refer to the following paper:

   M. Attene.
   A lightweight approach to repairing digitized polygon meshes.
   The Visual Computer, 2010. (c) Springer. DOI: 10.1007/s00371-010-0416-3

This software is based on ideas published therein. If you use MeshFix for research purposes you should cite the above paper in your published results. MeshFix cannot be used for commercial purposes without a proper licensing contract.

----------
PARAMETERS
----------
The user may want to force the software to join the connected components [-a] whose boundary loops are closer to each other. Otherwise, only the largest component is kept while the others are considered as "noise". By default, the output file is in OFF format. The user may change this to STL [-j]. If the software is run as part of a script that processes all the files in a directory, the parameter [-x] may be useful to avoid to re-run the repairing of a mesh if its result already exists as a file.

-------------------
COMMAND LINE SYNTAX
-------------------
Run MeshFix without parameters to get the synopsis.

---------------------
OTHER LAUNCHING MODES
---------------------
A mesh can also be fixed my simply dragging its icon on MeshFix's icon. In this case, only default parameters can be used.
MeshFix does not require any interaction, thus it can be inserted into a script to automatically repair all the models of a given repository (e.g. a folder). If some failure cases occur, they will be logged to "meshfix.log" and thus can be easily located and possibly processed is a second stage through interactive software tools such as ReMESH (http://remesh.sourceforge.net).

---------------------
SOURCE CODE
---------------------

From version 2.0 MeshFix is self-contained, that is, it does not depend on any other software.
To compile the source code on Windows you need Microsoft Visual C++ 2013 or newer.
Just click on vc12/MeshFix_All.sln and hit F7.
The executable will be saved in bin/.

Both 32bit and 64bit versions can be produced.

The source code is standard ANSI C++ and should be portable.
To compile on other operating systems (e.g. Linux) you may need to manually edit appropriate compilation instructions (e.g. Makefiles).


---------
Copyright
---------

MeshFix is

Copyright(C) 2010: IMATI-GE / CNR                                       

All rights reserved.                                                      
                                                                  
This program is dual-licensed as follows:

(1) You may use MeshFix as free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published 
by the Free Software Foundation; either version 3 of the License, or     
(at your option) any later version.                                      
In this case the program is distributed in the hope that it will be      
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)         
for more details.                                                        
                                                                         
(2) You may use MeshFix as part of a commercial software. In this case a
proper agreement must be reached with the Authors and with IMATI-GE/CNR  
based on a proper licensing contract.
