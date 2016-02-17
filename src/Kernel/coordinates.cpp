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

#include "coordinates.h"
#include <limits>
#include <cmath>


namespace T_MESH
{
#ifndef NAN
#define NAN std::numeric_limits<double>::quiet_NaN()
#endif

#ifdef USE_HYBRID_KERNEL

bool coord::use_rationals = false;

void coord::switchToDouble()
{
	if (_whv)
	{
		EXACT_NT * ov = (EXACT_NT *)_val;
		double d = (double)((EXACT_NT_DENOMINATOR(ov) != 0) ? (EXACT_NT_TO_DOUBLE((*ov))) : (NAN));
		_val = d2int64t(d);
		delete ov;
		_whv = false;
	}
}

void coord::switchToRational()
{
	if (!_whv)
	{
		double od = int64t2d(_val);
		if (od == NAN) _val = (int64_t)new EXACT_NT(0, 0);
		else _val = (int64_t)new EXACT_NT(od);
		_whv = true;
	}
}

 void coord::operator+=(const coord& a)
{
	if (use_rationals) { switchToRational(); getVal() += a.toRational(); }
	else { switchToDouble(); getDVal() += a.toDouble(); }
}

 void coord::operator-=(const coord& a)
{
	if (use_rationals) { switchToRational(); getVal() -= a.toRational(); }
	else { switchToDouble(); getDVal() -= a.toDouble(); }
}

 void coord::operator*=(const coord& a)
{
	if (use_rationals) { switchToRational(); getVal() *= a.toRational(); }
	else { switchToDouble(); getDVal() *= a.toDouble(); }
}

 void coord::operator/=(const coord& a)
{
	if (use_rationals) { switchToRational(); getVal() /= a.toRational(); }
	else { switchToDouble(); getDVal() /= a.toDouble(); }
}

 coord coord::operator+(const coord& a) const
{
	if (use_rationals) return coord(toRational() + a.toRational());
	else return coord(toDouble() + a.toDouble());
}

 coord coord::operator-(const coord& a) const
{
	if (use_rationals) return coord(toRational() - a.toRational());
	else return coord(toDouble() - a.toDouble());
}

 coord coord::operator*(const coord& a) const
{
	if (use_rationals) return coord(toRational() * a.toRational());
	else return coord(toDouble() * a.toDouble());
}

 coord coord::operator/(const coord& a) const
{
	if (use_rationals) return coord(toRational() / a.toRational());
	else return coord(toDouble() / a.toDouble());
}

 bool coord::operator==(const coord& a) const
{
	 if (_whv || a._whv /*use_rationals*/) return (toRational() == a.toRational());
	 else return (toDouble() == a.toDouble());
 }

 bool coord::operator!=(const coord& a) const
{
	 if (_whv || a._whv /*use_rationals*/) return (toRational() != a.toRational());
	 else return (toDouble() != a.toDouble());
}

 coord& coord::operator=(const coord& a)
 {
	 if (_whv) delete ((EXACT_NT *)_val);
	 _whv = a._whv; _val = (_whv) ? ((int64_t)new EXACT_NT(a.getVal())) : (a._val);
	 return *this;
 }

 void coord::setFromRational(const EXACT_NT& a)
 {
	 if (_whv) delete ((EXACT_NT *)_val);
	 _whv = 1; _val = (int64_t)new EXACT_NT(a);
 }

 bool coord::operator<(const coord& a) const
 {
	 if (_whv || a._whv /*use_rationals*/) return (toRational() < a.toRational());
	 else return (toDouble() < a.toDouble());
 }

bool coord::operator>(const coord& a) const
{
	if (_whv || a._whv /*use_rationals*/) return (toRational() > a.toRational());
	else return (toDouble() > a.toDouble());
}

 coord operator-(const coord& a)
{
	if (coord::use_rationals) return coord(-(a.toRational()));
	else return coord(-(a.toDouble()));
}

coord ceil(const coord& a)
{
	if (coord::use_rationals)
	{
		mpz_t n, d, f;
		mpz_init(n); mpz_init(d); mpz_init(f);
#ifdef USE_CGAL_LAZYNT
		mpz_set(n, a.toRational().exact().numerator().mpz());
		mpz_set(d, a.toRational().exact().denominator().mpz());
		mpz_cdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return coord(EXACT_NT(CGAL::Gmpz(f)));
#else
		mpz_set(n, a.toRational().get_num_mpz_t());
		mpz_set(d, a.toRational().get_den_mpz_t());
		mpz_cdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return coord(mpq_class(mpz_class(f)));
#endif
	}
	else
		return coord(::ceil(a.toDouble()));
}

coord floor(const coord& a)
{
	if (coord::use_rationals)
	{
		mpz_t n, d, f;
		mpz_init(n); mpz_init(d); mpz_init(f);
#ifdef USE_CGAL_LAZYNT
		mpz_set(n, a.toRational().exact().numerator().mpz());
		mpz_set(d, a.toRational().exact().denominator().mpz());
		mpz_fdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return coord(EXACT_NT(CGAL::Gmpz(f)));
#else
		mpz_set(n, a.toRational().get_num_mpz_t());
		mpz_set(d, a.toRational().get_den_mpz_t());
		mpz_fdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return coord(mpq_class(mpz_class(f)));
#endif
	}
	else
		return coord(::floor(a.toDouble()));
}

extern "C" double orient2d(double *, double *, double *);

coord coord::orient2D(const coord& px, const coord& py, const coord& qx, const coord& qy, const coord& rx, const coord& ry)
{
	if (use_rationals) return ((px - rx)*(qy - ry) - (py - ry)*(qx - rx));
	else
	{
		double pqr[6];
		pqr[0] = px.toDouble();  pqr[1] = py.toDouble();
		pqr[2] = qx.toDouble();  pqr[3] = qy.toDouble();
		pqr[4] = rx.toDouble();  pqr[5] = ry.toDouble();
		return orient2d(pqr, pqr + 2, pqr + 4);
	}
}

#endif
} //namespace T_MESH
