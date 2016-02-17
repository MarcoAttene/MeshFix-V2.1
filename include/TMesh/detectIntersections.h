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

#ifndef DETECT_INTERSECTIONS_H
#define DETECT_INTERSECTIONS_H

#include "tin.h"

namespace T_MESH
{
#define DI_MAX_NUMBER_OF_CELLS	10000
#define DI_EPSILON_POINT Point(1.0e-9, 1.0e-9, 1.0e-9)

class di_cell
{
public:
    Point mp, Mp;
	List triangles;

	di_cell() {}
	di_cell(Basic_TMesh *tin, bool useAll = true);

    bool is_triangleBB_in_cell(Triangle *t) const;

	di_cell *fork();
	void selectIntersections(bool justproper = false);
	bool doesNotIntersectForSure();
};

} //namespace T_MESH

#endif // DETECT_INTERSECTIONS_H
