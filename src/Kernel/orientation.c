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

// This code is inspired on ideas first published in the following paper:
// Jonathan Richard Shewchuk. Adaptive Precision Floating-Point Arithmetic 
// and Fast Robust Geometric Predicates, Discrete & Computational Geometry
// 18(3):305–363, October 1997.
//

#include <math.h>

#ifdef SPECIFY_FP_PRECISION
#include <float.h>
#endif

/*****************************************************************************/
/*                                                                           */
/*  This section contains private macros and functions.                      */
/*  These are not part of the public interface.                              */
/*                                                                           */
/*****************************************************************************/

#define FABS(a) (((a)>=0.0)?(a):(-(a)))
#define FTST(a,b,x,y) _bvr=x-a; y=b-_bvr
#define FTS(a,b,x,y) x=(double)(a+b); FTST(a,b,x,y)
#define TST(a,b,x,y) _bvr=(double)(x-a); _avr=x-_bvr; _brn=b-_bvr; _arn=a-_avr; y=_arn+_brn
#define TWS(a,b,x,y) x=(double)(a+b); TST(a,b,x,y)
#define TDT(a,b,x,y) _bvr=(double)(a-x); _avr=x+_bvr; _brn=_bvr-b; _arn=a-_avr; y=_arn+_brn
#define TWD(a,b,x,y) x=(double)(a-b); TDT(a,b,x,y)
#define SPLT(a,ahi,alo) c=(double)(_spl*a); abig=(double)(c-a); ahi=c-abig; alo=a-ahi
#define TPT(a,b,x,y) SPLT(a,ahi,alo); SPLT(b,bhi,blo); err1=x-(ahi*bhi); err2=err1-(alo*bhi); err3=err2-(ahi*blo); y=(alo*blo)-err3
#define TWP(a,b,x,y) x=(double)(a*b); TPT(a,b,x,y)
#define TPP(a,b,bhi,blo,x,y) x=(double)(a*b); SPLT(a,ahi,alo); err1=x-(ahi*bhi); err2=err1-(alo*bhi); err3=err2-(ahi*blo); y=(alo*blo)-err3
#define TOD(a1,a0,b,x2,x1,x0) TWD(a0,b,_i,x0); TWS(a1,_i,x2,x1)
#define TTD(a1,a0,b1,b0,x3,x2,x1,x0) TOD(a1,a0,b0,_j,_0,x0); TOD(_j,_0,b1,x3,x2,x1)
#define TOP(a1,a0,b,x3,x2,x1,x0) SPLT(b,bhi,blo); TPP(a0,b,bhi,blo,_i,x0); TPP(a1,b,bhi,blo,_j,_0); TWS(_i,_0,_k,x1); FTS(_j,_k,x3,x2)

double _spl, _eps, _reb, _ccwebA, _ccwebB, _ccwebC, _o3ebA, _o3ebB, _o3ebC;
double _iccebA, _iccebB, _iccebC, _ispebA, _ispebB, _ispebC;

int _fesze(int elen, double *e, int flen, double *f, double *h)
{
  double Q, Qnew, hh, _bvr, _avr, _brn, _arn, enow, fnow;
  int eindex, findex, hindex;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = e[++eindex];
  } else {
    Q = fnow;
    fnow = f[++findex];
  }
  hindex = 0;
  if ((eindex < elen) && (findex < flen)) {
    if ((fnow > enow) == (fnow > -enow)) {
      FTS(enow, Q, Qnew, hh);
      enow = e[++eindex];
    } else {
      FTS(fnow, Q, Qnew, hh);
      fnow = f[++findex];
    }
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
    while ((eindex < elen) && (findex < flen)) {
      if ((fnow > enow) == (fnow > -enow)) {
        TWS(Q, enow, Qnew, hh);
        enow = e[++eindex];
      } else {
        TWS(Q, fnow, Qnew, hh);
        fnow = f[++findex];
      }
      Q = Qnew;
      if (hh != 0.0) {
        h[hindex++] = hh;
      }
    }
  }
  while (eindex < elen) {
    TWS(Q, enow, Qnew, hh);
    enow = e[++eindex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  while (findex < flen) {
    TWS(Q, fnow, Qnew, hh);
    fnow = f[++findex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

int _seze(int elen, double *e, double b, double *h)
{
 double Q, sum, hh, product1, product0, enow, _bvr, _avr, _brn, _arn, c;
 double abig, ahi, alo, bhi, blo, err1, err2, err3;
  int eindex, hindex;

  SPLT(b, bhi, blo);
  TPP(e[0], b, bhi, blo, Q, hh);
  hindex = 0;
  if (hh != 0) {
    h[hindex++] = hh;
  }
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    TPP(enow, b, bhi, blo, product1, product0);
    TWS(Q, product0, sum, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
    FTS(product1, sum, Q, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

double _estm(int elen, double *e)
{
 int eindex;
 double Q = e[0];
 for (eindex = 1; eindex < elen; eindex++) Q += e[eindex];
 return Q;
}


double _adaptive2dorientation(double *pa, double *pb, double *pc, double detsum)
{
 double acx, acy, bcx, bcy,acxtail, acytail, bcxtail, bcytail, detleft, detright;
 double detlefttail, detrighttail, det, errbound, B[4], C1[8], C2[12], D[16];
 double B3, u[4], u3, s1, t1, s0, t0, _bvr, _avr, _brn, _arn, c;
 double abig, ahi, alo, bhi, blo, err1, err2, err3, _i, _j, _0;
 int C1length, C2length, Dlength;

  acx = (double) (pa[0] - pc[0]);
  bcx = (double) (pb[0] - pc[0]);
  acy = (double) (pa[1] - pc[1]);
  bcy = (double) (pb[1] - pc[1]);

  TWP(acx, bcy, detleft, detlefttail);
  TWP(acy, bcx, detright, detrighttail);

  TTD(detleft, detlefttail, detright, detrighttail,
               B3, B[2], B[1], B[0]);
  B[3] = B3;

  det = _estm(4, B);
  errbound = _ccwebB * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  TDT(pa[0], pc[0], acx, acxtail);
  TDT(pb[0], pc[0], bcx, bcxtail);
  TDT(pa[1], pc[1], acy, acytail);
  TDT(pb[1], pc[1], bcy, bcytail);

  if ((acxtail == 0.0) && (acytail == 0.0)
      && (bcxtail == 0.0) && (bcytail == 0.0)) {
    return det;
  }

  errbound = _ccwebC * detsum + _reb * FABS(det);
  det += (acx * bcytail + bcy * acxtail)
       - (acy * bcxtail + bcx * acytail);
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  TWP(acxtail, bcy, s1, s0);
  TWP(acytail, bcx, t1, t0);
  TTD(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C1length = _fesze(4, B, 4, u, C1);

  TWP(acx, bcytail, s1, s0);
  TWP(acy, bcxtail, t1, t0);
  TTD(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  C2length = _fesze(C1length, C1, 4, u, C2);

  TWP(acxtail, bcytail, s1, s0);
  TWP(acytail, bcxtail, t1, t0);
  TTD(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
  u[3] = u3;
  Dlength = _fesze(C2length, C2, 4, u, D);

  return(D[Dlength - 1]);
}

double _adaptive3dorientation(double *pa, double *pb, double *pc, double *pd, double permanent)
{
 double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz, det, errbound;
 double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
 double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0, bc[4], ca[4], ab[4];
 double bc3, ca3, ab3, adet[8], bdet[8], cdet[8];
 double abdet[16], *finnow, *finother, *finswap, fin1[192], fin2[192];
 double adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail, adztail, bdztail, cdztail;
 double at_blarge, at_clarge, bt_clarge, bt_alarge, ct_alarge, ct_blarge;
 double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
 double bdxt_cdy1, cdxt_bdy1, cdxt_ady1, adxt_cdy1, adxt_bdy1, bdxt_ady1;
 double bdxt_cdy0, cdxt_bdy0, cdxt_ady0, adxt_cdy0, adxt_bdy0, bdxt_ady0;
 double bdyt_cdx1, cdyt_bdx1, cdyt_adx1, adyt_cdx1, adyt_bdx1, bdyt_adx1;
 double bdyt_cdx0, cdyt_bdx0, cdyt_adx0, adyt_cdx0, adyt_bdx0, bdyt_adx0;
 double bct[8], cat[8], abt[8], bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
 double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1, bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
 double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0, u[4], v[12], w[16], u3, negate;
 double _bvr, _avr, _brn, _arn, c, abig, ahi, alo, bhi, blo;
 double err1, err2, err3, _i, _j, _k, _0;
 int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen, finlength;
 int vlength, wlength, alen, blen, clen, ablen, bctlen, catlen, abtlen;

  adx = (double) (pa[0] - pd[0]);
  bdx = (double) (pb[0] - pd[0]);
  cdx = (double) (pc[0] - pd[0]);
  ady = (double) (pa[1] - pd[1]);
  bdy = (double) (pb[1] - pd[1]);
  cdy = (double) (pc[1] - pd[1]);
  adz = (double) (pa[2] - pd[2]);
  bdz = (double) (pb[2] - pd[2]);
  cdz = (double) (pc[2] - pd[2]);

  TWP(bdx, cdy, bdxcdy1, bdxcdy0);
  TWP(cdx, bdy, cdxbdy1, cdxbdy0);
  TTD(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  alen = _seze(4, bc, adz, adet);

  TWP(cdx, ady, cdxady1, cdxady0);
  TWP(adx, cdy, adxcdy1, adxcdy0);
  TTD(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  blen = _seze(4, ca, bdz, bdet);

  TWP(adx, bdy, adxbdy1, adxbdy0);
  TWP(bdx, ady, bdxady1, bdxady0);
  TTD(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  clen = _seze(4, ab, cdz, cdet);

  ablen = _fesze(alen, adet, blen, bdet, abdet);
  finlength = _fesze(ablen, abdet, clen, cdet, fin1);

  det = _estm(finlength, fin1);
  errbound = _o3ebB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  TDT(pa[0], pd[0], adx, adxtail);
  TDT(pb[0], pd[0], bdx, bdxtail);
  TDT(pc[0], pd[0], cdx, cdxtail);
  TDT(pa[1], pd[1], ady, adytail);
  TDT(pb[1], pd[1], bdy, bdytail);
  TDT(pc[1], pd[1], cdy, cdytail);
  TDT(pa[2], pd[2], adz, adztail);
  TDT(pb[2], pd[2], bdz, bdztail);
  TDT(pc[2], pd[2], cdz, cdztail);

  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
      && (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
    return det;
  }

  errbound = _o3ebC * permanent + _reb * FABS(det);
  det += (adz * ((bdx * cdytail + cdy * bdxtail)
                 - (bdy * cdxtail + cdx * bdytail))
          + adztail * (bdx * cdy - bdy * cdx))
       + (bdz * ((cdx * adytail + ady * cdxtail)
                 - (cdy * adxtail + adx * cdytail))
          + bdztail * (cdx * ady - cdy * adx))
       + (cdz * ((adx * bdytail + bdy * adxtail)
                 - (ady * bdxtail + bdx * adytail))
          + cdztail * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if (adxtail == 0.0) {
    if (adytail == 0.0) {
      at_b[0] = 0.0;
      at_blen = 1;
      at_c[0] = 0.0;
      at_clen = 1;
    } else {
      negate = -adytail;
      TWP(negate, bdx, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      TWP(adytail, cdx, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    }
  } else {
    if (adytail == 0.0) {
      TWP(adxtail, bdy, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      negate = -adxtail;
      TWP(negate, cdy, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    } else {
      TWP(adxtail, bdy, adxt_bdy1, adxt_bdy0);
      TWP(adytail, bdx, adyt_bdx1, adyt_bdx0);
      TTD(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                   at_blarge, at_b[2], at_b[1], at_b[0]);
      at_b[3] = at_blarge;
      at_blen = 4;
      TWP(adytail, cdx, adyt_cdx1, adyt_cdx0);
      TWP(adxtail, cdy, adxt_cdy1, adxt_cdy0);
      TTD(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                   at_clarge, at_c[2], at_c[1], at_c[0]);
      at_c[3] = at_clarge;
      at_clen = 4;
    }
  }
  if (bdxtail == 0.0) {
    if (bdytail == 0.0) {
      bt_c[0] = 0.0;
      bt_clen = 1;
      bt_a[0] = 0.0;
      bt_alen = 1;
    } else {
      negate = -bdytail;
      TWP(negate, cdx, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      TWP(bdytail, adx, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    }
  } else {
    if (bdytail == 0.0) {
      TWP(bdxtail, cdy, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      negate = -bdxtail;
      TWP(negate, ady, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    } else {
      TWP(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
      TWP(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
      TTD(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                   bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
      bt_c[3] = bt_clarge;
      bt_clen = 4;
      TWP(bdytail, adx, bdyt_adx1, bdyt_adx0);
      TWP(bdxtail, ady, bdxt_ady1, bdxt_ady0);
      TTD(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                  bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
      bt_a[3] = bt_alarge;
      bt_alen = 4;
    }
  }
  if (cdxtail == 0.0) {
    if (cdytail == 0.0) {
      ct_a[0] = 0.0;
      ct_alen = 1;
      ct_b[0] = 0.0;
      ct_blen = 1;
    } else {
      negate = -cdytail;
      TWP(negate, adx, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      TWP(cdytail, bdx, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    }
  } else {
    if (cdytail == 0.0) {
      TWP(cdxtail, ady, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      negate = -cdxtail;
      TWP(negate, bdy, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    } else {
      TWP(cdxtail, ady, cdxt_ady1, cdxt_ady0);
      TWP(cdytail, adx, cdyt_adx1, cdyt_adx0);
      TTD(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                   ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
      ct_a[3] = ct_alarge;
      ct_alen = 4;
      TWP(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
      TWP(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
      TTD(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                   ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
      ct_b[3] = ct_blarge;
      ct_blen = 4;
    }
  }

  bctlen = _fesze(bt_clen, bt_c, ct_blen, ct_b, bct);
  wlength = _seze(bctlen, bct, adz, w);
  finlength = _fesze(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  catlen = _fesze(ct_alen, ct_a, at_clen, at_c, cat);
  wlength = _seze(catlen, cat, bdz, w);
  finlength = _fesze(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  abtlen = _fesze(at_blen, at_b, bt_alen, bt_a, abt);
  wlength = _seze(abtlen, abt, cdz, w);
  finlength = _fesze(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  if (adztail != 0.0) {
    vlength = _seze(4, bc, adztail, v);
    finlength = _fesze(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    vlength = _seze(4, ca, bdztail, v);
    finlength = _fesze(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    vlength = _seze(4, ab, cdztail, v);
    finlength = _fesze(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if (adxtail != 0.0) {
    if (bdytail != 0.0) {
      TWP(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
      TOP(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = _fesze(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        TOP(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = _fesze(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (cdytail != 0.0) {
      negate = -adxtail;
      TWP(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
      TOP(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = _fesze(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        TOP(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = _fesze(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (bdxtail != 0.0) {
    if (cdytail != 0.0) {
      TWP(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
      TOP(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = _fesze(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        TOP(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = _fesze(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (adytail != 0.0) {
      negate = -bdxtail;
      TWP(negate, adytail, bdxt_adyt1, bdxt_adyt0);
      TOP(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = _fesze(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        TOP(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = _fesze(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (cdxtail != 0.0) {
    if (adytail != 0.0) {
      TWP(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
      TOP(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = _fesze(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        TOP(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = _fesze(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (bdytail != 0.0) {
      negate = -cdxtail;
      TWP(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
      TOP(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = _fesze(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        TOP(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = _fesze(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }

  if (adztail != 0.0) {
    wlength = _seze(bctlen, bct, adztail, w);
    finlength = _fesze(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    wlength = _seze(catlen, cat, bdztail, w);
    finlength = _fesze(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    wlength = _seze(abtlen, abt, cdztail, w);
    finlength = _fesze(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  return finnow[finlength - 1];
}





/*****************************************************************************/
/*                                                                           */
/*  PUBLIC FUNCTIONS                                                         */
/*  initPredicates()   Sets  the variables used for exact arithmetic. This   */
/*                must be called once before using the other two functions.  */
/*  orient2d()    Computes the orientation of three 2D points.               */
/*  orient3d()    Computes the orientation of four 3D points.                */
/*                                                                           */
/*****************************************************************************/

void initPredicates()
{
 static char a_c=0;
 double hf, ck, lc;
 int e_o;

 if (a_c) return; else a_c = 1;

#ifdef SPECIFY_FP_PRECISION
 unsigned int old_cfp;
 _controlfp_s(&old_cfp, _PC_53, MCW_PC);
#endif
 
 e_o = 1;
 _eps = _spl = ck = 1.0;
 hf = 0.5;

 do
 {
  lc=ck;
  _eps *= hf;
  if (e_o) _spl *= 2.0;
  e_o = !e_o;
  ck = 1.0 + _eps;
 } while ((ck != 1.0) && (ck != lc));
 _spl += 1.0;

  _reb = (3.0 + 8.0 * _eps) * _eps;
  _ccwebA = (3.0 + 16.0 * _eps) * _eps;
  _ccwebB = (2.0 + 12.0 * _eps) * _eps;
  _ccwebC = (9.0 + 64.0 * _eps) * _eps * _eps;
  _o3ebA = (7.0 + 56.0 * _eps) * _eps;
  _o3ebB = (3.0 + 28.0 * _eps) * _eps;
  _o3ebC = (26.0 + 288.0 * _eps) * _eps * _eps;
  _iccebA = (10.0 + 96.0 * _eps) * _eps;
  _iccebB = (4.0 + 48.0 * _eps) * _eps;
  _iccebC = (44.0 + 576.0 * _eps) * _eps * _eps;
  _ispebA = (16.0 + 224.0 * _eps) * _eps;
  _ispebB = (5.0 + 72.0 * _eps) * _eps;
  _ispebC = (71.0 + 1408.0 * _eps) * _eps * _eps;

#ifdef SPECIFY_FP_PRECISION
  _controlfp_s(&old_cfp, _CW_DEFAULT, MCW_PC);
#endif
}

double orient2d(double *pa, double *pb, double *pc)
{
 double dlf, drg, det, dsm, eb;

 dlf = (pa[0]-pc[0])*(pb[1]-pc[1]);
 drg = (pa[1]-pc[1])*(pb[0]-pc[0]);
 det = dlf - drg;

 if (dlf > 0.0) {if (drg <= 0.0) return det; else dsm = dlf + drg;}
 else if (dlf < 0.0) {if (drg >= 0.0) return det; else dsm = -dlf - drg;}
 else return det;

 eb = _ccwebA*dsm;
 if ((det>=eb) || (-det>=eb)) return det;

 return _adaptive2dorientation(pa, pb, pc, dsm);
}

double orient3d(double *pa, double *pb, double *pc, double *pd)
{
 double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz, pm, eb;
 double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady, det;

 adx = pa[0]-pd[0]; bdx = pb[0]-pd[0]; cdx = pc[0]-pd[0];
 ady = pa[1]-pd[1]; bdy = pb[1]-pd[1]; cdy = pc[1]-pd[1];
 adz = pa[2]-pd[2]; bdz = pb[2]-pd[2]; cdz = pc[2]-pd[2];

 bdxcdy = bdx*cdy; cdxbdy = cdx*bdy;
 cdxady = cdx*ady; adxcdy = adx*cdy;
 adxbdy = adx*bdy; bdxady = bdx*ady;

 det = adz*(bdxcdy-cdxbdy)+bdz*(cdxady-adxcdy)+cdz*(adxbdy-bdxady);
 pm=(FABS(bdxcdy)+FABS(cdxbdy))*FABS(adz)+(FABS(cdxady)+FABS(adxcdy))*FABS(bdz)+(FABS(adxbdy)+FABS(bdxady))*FABS(cdz);
 eb = _o3ebA*pm;
 if ((det>eb) || (-det>eb)) return det;
 return _adaptive3dorientation(pa, pb, pc, pd, pm);
}
