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

#ifndef _COORDINATES_H
#define _COORDINATES_H

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#ifdef USE_HYBRID_KERNEL

#ifdef USE_CGAL_LAZYNT
#include <CGAL/Gmpq.h> 
#include <CGAL/Lazy_exact_nt.h> 
#include <CGAL/number_utils.h>

typedef CGAL::Lazy_exact_nt<CGAL::Gmpq> EXACT_NT;
#define EXACT_NT_TO_DOUBLE(x) (CGAL::to_double(x))
#define EXACT_NT_DENOMINATOR(x) ((x)->exact().denominator())
#define EXACT_NT_NUMERATOR(x) ((x)->exact().numerator())

#else
#include <mpirxx.h>

typedef mpq_class EXACT_NT;
#define EXACT_NT_TO_DOUBLE(x) ((x).get_d())
#define EXACT_NT_DENOMINATOR(x) ((x)->get_den())
#define EXACT_NT_NUMERATOR(x) ((x)->get_num())

#endif

#endif

#ifndef __APPLE__
#include <cstdint>
#endif

namespace T_MESH
{

#ifdef USE_HYBRID_KERNEL
class PM_Rational
	{
	private:
		// Only one of the following two static members can be true at any time, but both of them can be false
		// Thus, three configurations are possible:
		// 1) APPROXIMATED: use_rationals = false, use_filtering = false
		// 2) FILTERED: use_rationals = false, use_filtering = true
		// 3) PRECISE: use_rationals = true, use_filtering = false
		static bool use_rationals;	// Kernel_mode uses rational numbers (precise)
		static bool use_filtering;	// Kernel_mode uses filtering

	public:
		inline static bool isUsingRationals() { return use_rationals; }
		inline static bool isUsingFiltering() { return use_filtering; }
		static void useRationals(bool v) { use_rationals = v; use_filtering = !v; }
		static void useFiltering(bool v) { use_filtering = v; if (v) use_rationals = false; }

	protected:
		int64_t _val;		// value. Might contain either a double or a pointer to mpq_class.
		bool _whv;			// Which value type is stored here. 1=rational, 0=double

		inline bool isOfRationalType() const { return _whv; }
		inline bool isOfDoubleType() const { return !_whv; }

		inline static int64_t d2int64t(double a) { return *((int64_t *)((void *)(&a))); }
		inline static double& int64t2d(const int64_t& a) { return *((double *)((void *)(&a))); }

		inline EXACT_NT& getVal() { return *((EXACT_NT *)_val); }
		inline double& getDVal() { return int64t2d(_val); }

		inline const EXACT_NT& getVal() const { return *((EXACT_NT *)_val); }
		inline const double& getDVal() const { return int64t2d(_val); }

		void switchToDouble();
		void switchToRational();

	public:
		PM_Rational() : _whv(0) {} // Undetermined double
		PM_Rational(const EXACT_NT& a) { _whv = use_rationals; _val = (_whv) ? ((int64_t)new EXACT_NT(a)) : (d2int64t(EXACT_NT_TO_DOUBLE(a))); }
		PM_Rational(float a) { _whv = use_rationals; _val = (_whv) ? ((int64_t)new EXACT_NT(a)) : (d2int64t(a)); }
		PM_Rational(double a) { _whv = use_rationals; _val = (_whv) ? ((int64_t)new EXACT_NT(a)) : (d2int64t(a)); }
		PM_Rational(int a) { _whv = use_rationals; _val = (_whv) ? ((int64_t)new EXACT_NT(a)) : (d2int64t(a)); }

		PM_Rational(const PM_Rational& a) { _whv = a._whv; _val = (_whv) ? ((int64_t)new EXACT_NT(a.getVal())) : (a._val); }
		~PM_Rational() { if (_whv) delete ((EXACT_NT *)_val); }

		inline EXACT_NT toRational() const { return (_whv) ? (getVal()) : (EXACT_NT(getDVal())); }

		// The following toSomething() functions may cause rounding. Use with caution !
		inline double toDouble() const { return ((_whv) ? (EXACT_NT_TO_DOUBLE(getVal())) : (getDVal())); }
		inline int toInt() const { return int(toDouble()); }
		inline float toFloat() const { return float(toDouble()); }

		void operator+=(const PM_Rational& a);
		void operator-=(const PM_Rational& a);
		void operator*=(const PM_Rational& a);
		void operator/=(const PM_Rational& a);
		PM_Rational operator+(const PM_Rational& a) const;
		PM_Rational operator-(const PM_Rational& a) const;
		PM_Rational operator*(const PM_Rational& a) const;
		PM_Rational operator/(const PM_Rational& a) const;
		bool operator==(const PM_Rational& a) const;
		bool operator!=(const PM_Rational& a) const;

		PM_Rational& operator=(const PM_Rational& a);
		void setFromRational(const EXACT_NT& a);

		bool operator<(const PM_Rational& a) const;
		bool operator>(const PM_Rational& a) const;
		inline bool operator<=(const PM_Rational& a) const { return (operator==(a) || operator<(a)); }
		inline bool operator>=(const PM_Rational& a) const { return (operator==(a) || operator>(a)); }

		friend PM_Rational orient2D(const PM_Rational& px, const PM_Rational& py, const PM_Rational& qx, const PM_Rational& qy, const PM_Rational& rx, const PM_Rational& ry);
		friend PM_Rational orient3D(const class Point *t, const class Point *a, const class Point *b, const class Point *c);
};

PM_Rational operator-(const PM_Rational& a);
PM_Rational ceil(const PM_Rational& a);
PM_Rational floor(const PM_Rational& a);
PM_Rational round(const PM_Rational& a);

/**************** I/O operators ****************/

inline std::ostream & operator<<(std::ostream &o, const PM_Rational& c) { return o << c.toRational(); }
inline std::istream & operator>>(std::istream &i, PM_Rational& c) { EXACT_NT a; i >> a; c.setFromRational(a); return i; }

#define TMESH_TO_DOUBLE(x) ((x).toDouble())
#define TMESH_TO_FLOAT(x) ((x).toFloat())
#define TMESH_TO_INT(x) ((x).toInt())

#else

typedef double PM_Rational;

#define TMESH_TO_DOUBLE(x) (x)
#define TMESH_TO_FLOAT(x) ((float)(x))
#define TMESH_TO_INT(x) ((int)(x))

#endif

typedef PM_Rational coord;

#define TMESH_DETERMINANT3X3(a11, a12, a13, a21, a22, a23, a31, a32, a33) ((a11)*((a22)*(a33) - (a23)*(a32)) - (a12)*((a21)*(a33) - (a23)*(a31)) + (a13)*((a21)*(a32) - (a22)*(a31)))

} //namespace T_MESH

#endif //_COORDINATES_H

