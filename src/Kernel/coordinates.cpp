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

// Default behaviour = FILTERED KERNEL
bool PM_Rational::use_rationals = false;
bool PM_Rational::use_filtering = true;

void PM_Rational::switchToDouble()
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

void PM_Rational::switchToRational()
{
	if (!_whv)
	{
		double od = int64t2d(_val);
		if (od == NAN) _val = (int64_t)new EXACT_NT(0, 0);
		else _val = (int64_t)new EXACT_NT(od);
		_whv = true;
	}
}

 void PM_Rational::operator+=(const PM_Rational& a)
{
	if (use_rationals) { switchToRational(); getVal() += a.toRational(); }
	else { switchToDouble(); getDVal() += a.toDouble(); }
}

 void PM_Rational::operator-=(const PM_Rational& a)
{
	if (use_rationals) { switchToRational(); getVal() -= a.toRational(); }
	else { switchToDouble(); getDVal() -= a.toDouble(); }
}

 void PM_Rational::operator*=(const PM_Rational& a)
{
	if (use_rationals) { switchToRational(); getVal() *= a.toRational(); }
	else { switchToDouble(); getDVal() *= a.toDouble(); }
}

 void PM_Rational::operator/=(const PM_Rational& a)
{
	if (use_rationals) { switchToRational(); getVal() /= a.toRational(); }
	else { switchToDouble(); getDVal() /= a.toDouble(); }
}

 PM_Rational PM_Rational::operator+(const PM_Rational& a) const
{
	if (use_rationals) return PM_Rational(toRational() + a.toRational());
	else return PM_Rational(toDouble() + a.toDouble());
}

 PM_Rational PM_Rational::operator-(const PM_Rational& a) const
{
	if (use_rationals) return PM_Rational(toRational() - a.toRational());
	else return PM_Rational(toDouble() - a.toDouble());
}

 PM_Rational PM_Rational::operator*(const PM_Rational& a) const
{
	if (use_rationals) return PM_Rational(toRational() * a.toRational());
	else return PM_Rational(toDouble() * a.toDouble());
}

 PM_Rational PM_Rational::operator/(const PM_Rational& a) const
{
	if (use_rationals) return PM_Rational(toRational() / a.toRational());
	else return PM_Rational(toDouble() / a.toDouble());
}

 bool PM_Rational::operator==(const PM_Rational& a) const
{
	 if (_whv || a._whv /*use_rationals*/) return (toRational() == a.toRational());
	 else return (toDouble() == a.toDouble());
 }

 bool PM_Rational::operator!=(const PM_Rational& a) const
{
	 if (_whv || a._whv /*use_rationals*/) return (toRational() != a.toRational());
	 else return (toDouble() != a.toDouble());
}

 PM_Rational& PM_Rational::operator=(const PM_Rational& a)
 {
	 if (_whv) delete ((EXACT_NT *)_val);
	 _whv = a._whv; _val = (_whv) ? ((int64_t)new EXACT_NT(a.getVal())) : (a._val);
	 return *this;
 }

 void PM_Rational::setFromRational(const EXACT_NT& a)
 {
	 if (_whv) delete ((EXACT_NT *)_val);
	 _whv = 1; _val = (int64_t)new EXACT_NT(a);
 }

 bool PM_Rational::operator<(const PM_Rational& a) const
 {
	 if (_whv || a._whv) return (toRational() < a.toRational());
	 else return (toDouble() < a.toDouble());
 }

bool PM_Rational::operator>(const PM_Rational& a) const
{
	if (_whv || a._whv) return (toRational() > a.toRational());
	else return (toDouble() > a.toDouble());
}

PM_Rational operator-(const PM_Rational& a) // This might be probably changed... do not understand why to switch..
{
	if (PM_Rational::isUsingRationals()) return PM_Rational(-(a.toRational()));
	else return PM_Rational(-(a.toDouble()));
}

PM_Rational ceil(const PM_Rational& a)
{
	if (PM_Rational::isUsingRationals())
	{
		mpz_t n, d, f;
		mpz_init(n); mpz_init(d); mpz_init(f);
#ifdef USE_CGAL_LAZYNT
		mpz_set(n, a.toRational().exact().numerator().mpz());
		mpz_set(d, a.toRational().exact().denominator().mpz());
		mpz_cdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return PM_Rational(EXACT_NT(CGAL::Gmpz(f)));
#else
		mpz_set(n, a.toRational().get_num_mpz_t());
		mpz_set(d, a.toRational().get_den_mpz_t());
		mpz_cdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return PM_Rational(mpq_class(mpz_class(f)));
#endif
	}
	else
		return PM_Rational(::ceil(a.toDouble()));
}

PM_Rational floor(const PM_Rational& a)
{
	if (PM_Rational::isUsingRationals())
	{
		mpz_t n, d, f;
		mpz_init(n); mpz_init(d); mpz_init(f);
#ifdef USE_CGAL_LAZYNT
		mpz_set(n, a.toRational().exact().numerator().mpz());
		mpz_set(d, a.toRational().exact().denominator().mpz());
		mpz_fdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return PM_Rational(EXACT_NT(CGAL::Gmpz(f)));
#else
		mpz_set(n, a.toRational().get_num_mpz_t());
		mpz_set(d, a.toRational().get_den_mpz_t());
		mpz_fdiv_q(f, n, d);
		mpz_clear(n); mpz_clear(d);
		return PM_Rational(mpq_class(mpz_class(f)));
#endif
	} else
		return PM_Rational(::floor(a.toDouble()));
}

PM_Rational round(const PM_Rational& a)
{
	if (PM_Rational::isUsingRationals())
	{
		mpz_t n, d, f, c;
		mpz_init(n); mpz_init(d); mpz_init(f); mpz_init(c);
#ifdef USE_CGAL_LAZYNT
		mpz_set(n, a.toRational().exact().numerator().mpz());
		mpz_set(d, a.toRational().exact().denominator().mpz());
		mpz_fdiv_q(f, n, d);
		mpz_cdiv_q(c, n, d);
		mpz_clear(n); mpz_clear(d);
		PM_Rational fr = PM_Rational(EXACT_NT(CGAL::Gmpz(f)));
		PM_Rational cr = PM_Rational(EXACT_NT(CGAL::Gmpz(c)));
		mpz_clear(f); mpz_clear(c);
		return ((a - fr) < (cr - a)) ? (fr) : (cr);
#else
		mpz_set(n, a.toRational().get_num_mpz_t());
		mpz_set(d, a.toRational().get_den_mpz_t());
		mpz_fdiv_q(f, n, d);
		mpz_cdiv_q(c, n, d);
		mpz_clear(n); mpz_clear(d);
		PM_Rational fr = PM_Rational(mpq_class(mpz_class(f)));
		PM_Rational cr = PM_Rational(mpq_class(mpz_class(c)));
		mpz_clear(f); mpz_clear(c);
		return ((a - fr) < (cr - a)) ? (fr) : (cr);
#endif
	} else
		return PM_Rational(::round(a.toDouble()));
}

#endif
} //namespace T_MESH
