#include "slicot.h"
#include "slicot_blas.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void ib01pd(const char *meth, const char *job, const char *jobcv,
            i32 nobr, i32 n, i32 m, i32 l, i32 nsmpl,
            f64 *r, i32 ldr, f64 *a, i32 lda, f64 *c, i32 ldc,
            f64 *b, i32 ldb, f64 *d, i32 ldd, f64 *q, i32 ldq,
            f64 *ry, i32 ldry, f64 *s, i32 lds, f64 *o, i32 ldo,
            f64 tol, i32 *iwork, f64 *dwork, i32 ldwork,
            i32 *iwarn, i32 *info)
{
    const f64 ZERO = 0.0;
    const f64 ONE = 1.0;
    const f64 TWO = 2.0;
    const f64 THREE = 3.0;

    bool moesp = (meth[0] == 'M' || meth[0] == 'm');
    bool n4sid = (meth[0] == 'N' || meth[0] == 'n');
    bool withal = (job[0] == 'A' || job[0] == 'a');
    bool withc = (job[0] == 'C' || job[0] == 'c') || withal;
    bool withd = (job[0] == 'D' || job[0] == 'd') || withal;
    bool withb = (job[0] == 'B' || job[0] == 'b') || withd;
    bool withco = (jobcv[0] == 'C' || jobcv[0] == 'c');

    i32 mnobr = m * nobr;
    i32 lnobr = l * nobr;
    i32 lmnobr = lnobr + mnobr;
    i32 lmmnob = lnobr + 2 * mnobr;
    i32 mnobrn = mnobr + n;
    i32 lnobrn = lnobr - n;
    i32 ldun2 = lnobr - l;
    i32 ldunn = ldun2 * n;
    i32 lmmnol = lmmnob + l;
    i32 nr = lmnobr + lmnobr;
    i32 npl = n + l;
    i32 n2 = n + n;
    i32 nn = n * n;

    *iwarn = 0;
    *info = 0;

    if (!(moesp || n4sid)) {
        *info = -1;
    } else if (!(withb || withc)) {
        *info = -2;
    } else if (!(withco || jobcv[0] == 'N' || jobcv[0] == 'n')) {
        *info = -3;
    } else if (nobr <= 1) {
        *info = -4;
    } else if (n <= 0 || n >= nobr) {
        *info = -5;
    } else if (m < 0) {
        *info = -6;
    } else if (l <= 0) {
        *info = -7;
    } else if (withco && nsmpl < nr) {
        *info = -8;
    } else if (ldr < nr) {
        *info = -10;
    } else if (lda < 1 || ((withc || (withb && n4sid)) && lda < n)) {
        *info = -12;
    } else if (ldc < 1 || ((withc || (withb && n4sid)) && ldc < l)) {
        *info = -14;
    } else if (ldb < 1 || (withb && ldb < n && m > 0)) {
        *info = -16;
    } else if (ldd < 1 || (withd && ldd < l && m > 0)) {
        *info = -18;
    } else if (ldq < 1 || (withco && ldq < n)) {
        *info = -20;
    } else if (ldry < 1 || (withco && ldry < l)) {
        *info = -22;
    } else if (lds < 1 || (withco && lds < n)) {
        *info = -24;
    } else if (ldo < 1 || ((withco || n4sid) && ldo < lnobr)) {
        *info = -26;
    } else {
        i32 iaw = 0;
        i32 minwrk = ldunn + 4 * n;

        if (moesp) {
            if (withc) {
                i32 alt1 = 2 * ldunn + n2;
                i32 alt2 = ldunn + nn + 7 * n;
                if (alt1 > minwrk) minwrk = alt1;
                if (alt2 > minwrk) minwrk = alt2;
            }
        }

        if ((m > 0 && withb) || n4sid) {
            i32 id = n4sid ? n : 0;
            i32 alt = 2 * ldunn + nn + id + 7 * n;
            if (alt > minwrk) minwrk = alt;
            if (moesp) {
                i32 alt1 = ldunn + n + 6 * mnobr;
                i32 alt2 = ldunn + n + (l + mnobr > lnobr + (3 * lnobr + 1 > m ? 3 * lnobr + 1 : m)
                          ? l + mnobr : lnobr + (3 * lnobr + 1 > m ? 3 * lnobr + 1 : m));
                if (alt1 > minwrk) minwrk = alt1;
                if (alt2 > minwrk) minwrk = alt2;
            }
        } else {
            if (moesp)
                iaw = n + nn;
        }

        if (n4sid || withco) {
            i32 alt1 = ldunn + iaw + n2 + (5 * n > lmmnol ? 5 * n : lmmnol);
            i32 id = n4sid ? n : 0;
            i32 alt2 = id + 4 * mnobrn + 1;
            i32 alt3 = id + mnobrn + npl;
            if (alt1 > minwrk) minwrk = alt1;
            if (alt2 > minwrk) minwrk = alt2;
            if (alt3 > minwrk) minwrk = alt3;

            if (n4sid && (m > 0 && withb)) {
                i32 alt = mnobr * npl * (m * npl + 1) + (npl * npl > 4 * m * npl + 1 ? npl * npl : 4 * m * npl + 1);
                if (alt > minwrk) minwrk = alt;
            }
        }

        if (ldwork < minwrk) {
            *info = -30;
            dwork[0] = (f64)minwrk;
        }
    }

    if (*info != 0) {
        return;
    }

    i32 nr2 = mnobr;
    i32 nr3 = lmnobr;
    i32 nr4 = lmmnob;

    f64 eps = SLC_DLAMCH("Precision");
    f64 thresh = pow(eps, TWO / THREE);
    f64 svlmax = ZERO;
    f64 rcond4 = ONE;
    f64 rcond1, rcond2, rcond3;

    i32 igal = 0;
    SLC_DLACPY("Full", &ldun2, &n, &r[nr2 + nr2 * ldr], &ldr, &dwork[igal], &ldun2);

    i32 itau1 = igal + ldunn;
    i32 jwork = itau1 + n;
    i32 ldw = jwork;

    SLC_DGEQRF(&ldun2, &n, &dwork[igal], &ldun2, &dwork[itau1],
               &dwork[jwork], &(i32){ldwork - jwork}, info);

    SLC_DTRCON("1-norm", "Upper", "NonUnit", &n, &dwork[igal], &ldun2,
               &rcond1, &dwork[jwork], iwork, info);

    f64 toll1 = tol;
    if (toll1 <= ZERO) toll1 = (f64)(nn) * eps;

    if (moesp) {
        i32 ncol = 0;
        i32 iun2 = jwork;

        if (withc) {
            SLC_DLACPY("Full", &l, &n, &r[nr2 + nr2 * ldr], &ldr, c, &ldc);
            SLC_DLACPY("Full", &ldun2, &n, &r[nr2 + l + nr2 * ldr], &ldr, &dwork[iun2], &ldun2);

            jwork = iun2 + ldunn;
            SLC_DORMQR("Left", "Transpose", &ldun2, &n, &n, &dwork[igal], &ldun2,
                       &dwork[itau1], &dwork[iun2], &ldun2, &dwork[jwork],
                       &(i32){ldwork - jwork}, info);
            SLC_DLACPY("Full", &n, &n, &dwork[iun2], &ldun2, a, &lda);
            ncol = n;
            jwork = iun2;
        }

        i32 rank = n;
        if (rcond1 > (toll1 > thresh ? toll1 : thresh)) {
            if (withc) {
                SLC_DTRTRS("Upper", "NoTranspose", "NonUnit", &n, &n,
                           &dwork[igal], &ldun2, a, &lda, info);
                if (*info > 0) {
                    *info = 3;
                    return;
                }
            }
            rank = n;
        } else {
            i32 iu = iun2;
            i32 isv = iu + nn;
            jwork = isv + n;

            i32 ihous = igal;
            if (m > 0 && withb) {
                ihous = jwork;
                jwork = ihous + ldunn;
                SLC_DLACPY("Lower", &ldun2, &n, &dwork[igal], &ldun2, &dwork[ihous], &ldun2);
            }

            char jobp = (m > 0 && withb) ? 'P' : 'N';
            mb02ud("Not factored", "Left", "NoTranspose", &jobp, n, ncol, ONE, toll1,
                   &rank, &dwork[igal], ldun2, &dwork[iu], n, &dwork[isv],
                   a, lda, &r[nr3 + nr2 * ldr], ldr, &dwork[jwork], ldwork - jwork, info);

            if (*info != 0) {
                *info = 2;
                return;
            }

            if (rank == 0) {
                jobp = 'N';
            } else if (m > 0 && withb) {
                i32 diff = ldun2 - n;
                SLC_DLASET("Full", &n, &diff, &ZERO, &ZERO, &r[nr3 + (nr2 + n) * ldr], &ldr);
                SLC_DORMQR("Right", "Transpose", &n, &ldun2, &n, &dwork[ihous], &ldun2,
                           &dwork[itau1], &r[nr3 + nr2 * ldr], &ldr, &dwork[jwork],
                           &(i32){ldwork - jwork}, info);
                jwork = iun2;
            }
            ldw = jwork;
        }

        rcond2 = ONE;

        if (m > 0 && withb) {
            SLC_DTRCON("1-norm", "Upper", "NonUnit", &mnobr, &r[nr3], &ldr,
                       &rcond2, &dwork[jwork], iwork, info);

            f64 toll = tol;
            if (toll <= ZERO) toll = (f64)(mnobr * mnobr) * eps;

            SLC_DGEMM("Transpose", "Transpose", &lnobrn, &mnobr, &lnobr,
                      &ONE, &r[nr2 + (nr2 + n) * ldr], &ldr, &r[nr2 * ldr], &ldr,
                      &ZERO, &r[nr2 + nr3 * ldr], &ldr);

            if (rcond2 > (toll > thresh ? toll : thresh)) {
                SLC_DTRSM("Right", "Upper", "Transpose", "Non-unit", &lnobrn, &mnobr,
                          &ONE, &r[nr3], &ldr, &r[nr2 + nr3 * ldr], &ldr);
            } else {
                i32 isv = ldw;
                jwork = isv + mnobr;
                i32 rank11;
                mb02ud("Not factored", "Right", "Transpose", "No pinv", lnobrn, mnobr,
                       ONE, toll, &rank11, &r[nr3], ldr, &r[nr3 + nr3 * ldr], ldr,
                       &dwork[isv], &r[nr2 + nr3 * ldr], ldr, NULL, 1, &dwork[jwork],
                       ldwork - jwork, info);
                if (*info != 0) {
                    *info = 2;
                    return;
                }
                jwork = ldw;
            }

            char jobpy = withal ? 'D' : job[0];
            if (withco) {
                SLC_DLACPY("Full", &lnobr, &n, &r[nr2 + nr2 * ldr], &ldr, o, &ldo);
            }

            ib01py(meth, &jobpy, nobr, n, m, l, rank, &r[nr2 + nr2 * ldr], ldr,
                   &dwork[igal], ldun2, &dwork[itau1], &r[nr3 + nr2 * ldr], ldr,
                   &r[nr2 + nr3 * ldr], ldr, &r[nr4 + nr2 * ldr], ldr,
                   &r[nr4 + nr3 * ldr], ldr, b, ldb, d, ldd, tol, iwork,
                   &dwork[jwork], ldwork - jwork, iwarn, info);

            if (*info != 0) return;
            rcond4 = dwork[jwork + 1];

            if (withco) {
                SLC_DLACPY("Full", &lnobr, &n, o, &ldo, &r[nr2], &ldr);
            }
        } else {
            rcond2 = ONE;
        }

        if (!withco) {
            rcond3 = ONE;
            goto finish;
        }
    } else {
        rcond2 = ONE;
    }

    if (n4sid || (!moesp && withco)) {
        if (n4sid) {
            SLC_DLACPY("Full", &lnobr, &n, &r[nr2 + nr2 * ldr], &ldr, o, &ldo);
        }
    }

    rcond3 = ONE;

finish:
    dwork[0] = (f64)(ldwork);
    dwork[1] = rcond1;
    dwork[2] = rcond2;
    dwork[3] = rcond3;
    dwork[4] = rcond4;
}
