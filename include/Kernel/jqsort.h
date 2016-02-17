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

//! \file
//! \brief Declaration of a generic QuickSort function.
//!
//! The  jqsort()  function sorts an array with numels elements.
//! The v argument points to the start of the array of elements casted to void *.
//! The contents of the array are sorted in ascending order according to  a
//! comparison  function  pointed  to  by  comp, which is called with two
//! arguments that point to the objects being compared.
//! The comparison function must return an integer less than, equal to,  or
//! greater  than  zero  if  the first argument is considered to be respectively
//! less than, equal to, or greater than the second.  If two members
//! compare as equal, their order in the sorted array is undefined.
//! See the manpage of the standard library qsort() function for further information.

namespace T_MESH
{

extern void jqsort(void *v[], int numels, int (*comp)(const void *, const void *));

} //namespace T_MESH
