#include "slicot.h"
#include "slicot_blas.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

void sg03bv(const char* trans, const i32 n, const f64* a, const i32 lda,
            const f64* e, const i32 lde, f64* b, const i32 ldb,
            f64* scale, f64* dwork, i32* info)
{
    const f64 mone = -1.0;
    const f64 one = 1.0;
    const f64 two = 2.0;
    const f64 zero = 0.0;

    f64 m1[4], m2[4], tm[4], ui[4];

    f64 bignum, c, delta1, eps, r, s, scale1, smlnum, sqtwo, t, x;
    i32 i, info1, j, kb, kh, kl, kl1, l, ldws, uiipt, wpt, ypt;
    bool notrns;

    const i32 int0 = 0;
    const i32 int1 = 1;
    const i32 int2 = 2;

    notrns = (*trans == 'N' || *trans == 'n');

    *info = 0;
    if (!notrns && *trans != 'T' && *trans != 't') {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < (n > 0 ? n : 1)) {
        *info = -4;
    } else if (lde < (n > 0 ? n : 1)) {
        *info = -6;
    } else if (ldb < (n > 0 ? n : 1)) {
        *info = -8;
    }

    if (*info != 0) {
        return;
    }

    *scale = one;

    if (n == 0) {
        return;
    }

    sqtwo = sqrt(two);
    eps = SLC_DLAMCH("P");
    smlnum = SLC_DLAMCH("S") / eps;
    bignum = one / smlnum;
    SLC_DLABAD(&smlnum, &bignum);

    printf("\n=== SG03BV C DEBUG ===\n");
    printf("INPUT: N=%d, TRANS=%c\n", n, *trans);
    printf("INPUT B (first 3x3):\n");
    for (i32 ii = 0; ii < (n < 3 ? n : 3); ii++) {
        for (i32 jj = 0; jj < (n < 3 ? n : 3); jj++) {
            printf(" %23.16E", b[ii + jj*ldb]);
        }
        printf("\n");
    }
    printf("sqtwo = %.16E\n", sqtwo);
    printf("eps = %.16E\n", eps);
    printf("smlnum = %.16E\n", smlnum);
    printf("bignum = %.16E\n", bignum);

    uiipt = 0;
    wpt = 2 * n - 2;
    ypt = 4 * n - 4;
    ldws = n - 1;

    if (notrns) {
        kh = 0;

        while (kh < n) {
            kl = kh;
            if (kl == n - 1) {
                kh = n - 1;
                kb = 1;
            } else {
                if (kl + 1 < n && a[(kl + 1) + kl * lda] == zero) {
                    kh = kl;
                    kb = 1;
                } else {
                    kh = kl + 1;
                    kb = 2;
                }
            }

            printf("\n========== ITERATION KL=%d, KH=%d, KB=%d ==========\n", kl, kh, kb);

            printf("\nDEBUG: B diagonal before KB block processing:\n");
            for (i32 ii = 0; ii < n; ii++) {
                printf("  B[%d,%d] = %.16E\n", ii, ii, b[ii + ii*ldb]);
            }

            if (kb == 1) {
                printf("DEBUG: KB=1 block at KL=%d\n", kl);
                printf("  a[%d,%d] = %.16E\n", kl, kl, a[kl + kl*lda]);
                printf("  e[%d,%d] = %.16E\n", kl, kl, e[kl + kl*lde]);
                printf("  b[%d,%d] = %.16E\n", kl, kl, b[kl + kl*ldb]);

                delta1 = -a[kl * lda + kl];
                printf("  delta1 (before sqrt) = %.16E\n", delta1);
                if (delta1 <= zero) {
                    *info = 3;
                    return;
                }
                delta1 = sqrt(delta1) * sqrt(e[kl * lde + kl]);
                printf("  delta1 (after sqrt) = %.16E\n", delta1);
                t = (b[kl * ldb + kl] * smlnum) / sqtwo;
                printf("  t = %.16E\n", t);
                if (t > delta1) {
                    scale1 = delta1 / t;
                    printf("  SCALING: scale1 = %.16E\n", scale1);
                    *scale = scale1 * (*scale);
                    SLC_DLASCL("Upper", &int0, &int0, &one, &scale1, &n, &n,
                              b, &ldb, &info1);
                }
                ui[0] = b[kl * ldb + kl] / delta1 / sqtwo;
                printf("  ui[0] = %.16E\n", ui[0]);
                m1[0] = a[kl * lda + kl] / e[kl * lde + kl];
                printf("  m1[0] = %.16E\n", m1[0]);
                m2[0] = (delta1 / e[kl * lde + kl]) * sqtwo;
                printf("  m2[0] = %.16E\n", m2[0]);

            } else {
                printf("DEBUG: KB=2 block at KL=%d, KH=%d\n", kl, kh);
                printf("  Before SG03BX:\n");
                printf("    a[%d,%d] = %.16E, a[%d,%d] = %.16E\n", kl, kl, a[kl + kl*lda], kl, kl+1, a[kl + (kl+1)*lda]);
                printf("    a[%d,%d] = %.16E, a[%d,%d] = %.16E\n", kl+1, kl, a[(kl+1) + kl*lda], kl+1, kl+1, a[(kl+1) + (kl+1)*lda]);
                printf("    b[%d,%d] = %.16E, b[%d,%d] = %.16E\n", kl, kl, b[kl + kl*ldb], kl, kl+1, b[kl + (kl+1)*ldb]);
                printf("    b[%d,%d] = %.16E, b[%d,%d] = %.16E\n", kl+1, kl, b[(kl+1) + kl*ldb], kl+1, kl+1, b[(kl+1) + (kl+1)*ldb]);

                sg03bx("C", "N", &a[kl * lda + kl], lda, &e[kl * lde + kl],
                      lde, &b[kl * ldb + kl], ldb, ui, int2, &scale1,
                      m1, int2, m2, int2, info);

                if (*info != 0) {
                    return;
                }

                printf("  After SG03BX: scale1 = %.16E\n", scale1);
                printf("    ui = [%.16E, %.16E; %.16E, %.16E]\n", ui[0], ui[2], ui[1], ui[3]);
                printf("    m1 = [%.16E, %.16E; %.16E, %.16E]\n", m1[0], m1[2], m1[1], m1[3]);
                printf("    m2 = [%.16E, %.16E; %.16E, %.16E]\n", m2[0], m2[2], m2[1], m2[3]);

                if (scale1 != one) {
                    *scale = scale1 * (*scale);
                    SLC_DLASCL("Upper", &int0, &int0, &one, &scale1, &n, &n,
                              b, &ldb, &info1);
                    printf("  B scaled by %.16E\n", scale1);
                }
            }

            if (kh < n - 1) {
                i32 nkh = n - kh - 1;
                printf("\nDEBUG: Computing Y matrix (kh=%d < n-1=%d, nkh=%d)\n", kh, n-1, nkh);

                SLC_DGEMM("T", "N", &nkh, &kb, &kb, &mone, &b[kl + (kh + 1) * ldb],
                         &ldb, m2, &int2, &zero, &dwork[uiipt], &ldws);
                printf("  After B*m2 -> uiipt\n");

                SLC_DGEMM("T", "T", &nkh, &kb, &kb, &mone, &a[kl + (kh + 1) * lda],
                         &lda, ui, &int2, &one, &dwork[uiipt], &ldws);
                printf("  After A*ui added to uiipt\n");

                SLC_DGEMM("T", "N", &kb, &kb, &kb, &one, ui, &int2, m1, &int2,
                         &zero, tm, &int2);
                printf("  After ui*m1 -> tm\n");

                SLC_DGEMM("T", "N", &nkh, &kb, &kb, &mone, &e[kl + (kh + 1) * lde],
                         &lde, tm, &int2, &one, &dwork[uiipt], &ldws);
                printf("  After E*tm added to uiipt (Y computed)\n");

                SLC_DLASET("A", &kb, &kb, &zero, &one, tm, &int2);

                sg03bw("N", nkh, kb, &a[(kh + 1) * lda + kh + 1], lda,
                      tm, int2, &e[(kh + 1) * lde + kh + 1], lde, m1, int2,
                      &dwork[uiipt], ldws, &scale1, info);

                if (scale1 != one) {
                    *scale = scale1 * (*scale);
                    SLC_DLASCL("Upper", &int0, &int0, &one, &scale1, &n, &n,
                              b, &ldb, &info1);
                    i32 int4 = 4;
                    SLC_DSCAL(&int4, &scale1, ui, &int1);
                }

                SLC_DLACPY("A", &nkh, &kb, &dwork[uiipt], &ldws,
                          &dwork[wpt], &ldws);

                SLC_DTRMM("L", "U", "T", "N", &nkh, &kb, &one,
                         &e[(kh + 1) * lde + kh + 1], &lde, &dwork[wpt], &ldws);

                SLC_DGEMM("T", "T", &nkh, &kb, &kb, &one, &e[kl + (kh + 1) * lde],
                         &lde, ui, &int2, &one, &dwork[wpt], &ldws);

                SLC_DCOPY(&nkh, &b[(kh + 1) * ldb + kl], &ldb, &dwork[ypt], &int1);
                if (kh > kl) {
                    SLC_DCOPY(&nkh, &b[(kh + 1) * ldb + kh], &ldb,
                             &dwork[ypt + ldws], &int1);
                }

                printf("\nDEBUG: Y vector (dwork[ypt]) before wpt*m2:\n");
                for (i32 ii = 0; ii < nkh; ii++) {
                    printf("  Y[%d] = %.16E\n", ii, dwork[ypt + ii]);
                }

                printf("\nDEBUG: WPT matrix (nkh=%d, kb=%d):\n", nkh, kb);
                for (i32 ii = 0; ii < nkh; ii++) {
                    printf("  ");
                    for (i32 jj = 0; jj < kb; jj++) {
                        printf(" %.16E", dwork[wpt + ii + jj*ldws]);
                    }
                    printf("\n");
                }

                printf("\nDEBUG: M2 matrix (kb=%d):\n", kb);
                for (i32 ii = 0; ii < kb; ii++) {
                    printf("  ");
                    for (i32 jj = 0; jj < kb; jj++) {
                        printf(" %.16E", m2[ii + jj*2]);
                    }
                    printf("\n");
                }

                SLC_DGEMM("N", "T", &nkh, &kb, &kb, &mone, &dwork[wpt],
                         &ldws, m2, &int2, &one, &dwork[ypt], &ldws);

                printf("\nDEBUG: Y vector (dwork[ypt]) after wpt*m2 (final Y):\n");
                for (i32 ii = 0; ii < nkh; ii++) {
                    printf("  Y[%d] = %.16E\n", ii, dwork[ypt + ii]);
                }
                printf("DEBUG: ypt offset=%d, ldws=%d\n", ypt, ldws);

                printf("\nDEBUG QR: Starting QR factorization (nkh=%d, kb=%d)\n", nkh, kb);
                printf("  B diagonal before QR:\n");
                for (i32 ii = kh+1; ii < n; ii++) {
                    printf("    B[%d,%d] = %.16E\n", ii, ii, b[ii + ii*ldb]);
                }
                l = ypt;
                for (j = 0; j < kb; j++) {
                    printf("\n--- QR ITERATION J=%d (processing column %d) ---\n", j, j);
                    for (i = 0; i < nkh; i++) {
                        x = b[(kh + i + 1) + (kh + i + 1) * ldb];
                        t = dwork[l + i];
                        printf("  I=%d: X=b[%d,%d]=%.16E, T=dwork[%d]=%.16E\n",
                               i, kh+i+1, kh+i+1, x, l+i, t);
                        SLC_DLARTG(&x, &t, &c, &s, &r);
                        printf("       DLARTG output: R=%.16E, C=%.16E, S=%.16E\n", r, c, s);
                        b[(kh + i + 1) + (kh + i + 1) * ldb] = r;
                        printf("       Set b[%d,%d] = %.16E\n", kh+i+1, kh+i+1, r);
                        if (i < nkh - 1) {
                            i32 nmihm1 = nkh - i - 1;
                            printf("       Applying DROT to %d elements (b[%d:%d,%d] and dwork[%d:%d])\n",
                                   nmihm1, kh+i+2, kh+nkh, kh+i+1, l+i+1, l+nkh-1);
                            SLC_DROT(&nmihm1, &b[(kh + i + 2) * ldb + kh + i + 1],
                                    &ldb, &dwork[l + i + 1], &int1, &c, &s);
                            printf("       After DROT, b[%d,%d]=%.16E\n", kh+i+2, kh+i+1,
                                   (i < nkh-1) ? b[(kh+i+2) + (kh+i+1)*ldb] : 0.0);
                        }
                    }
                    printf("  After J=%d iteration, B diagonal:\n", j);
                    for (i32 ii = kh+1; ii < n; ii++) {
                        printf("    B[%d,%d] = %.16E\n", ii, ii, b[ii + ii*ldb]);
                    }
                    l += ldws;
                }

                for (i = kh + 1; i < n; i++) {
                    if (b[i + i * ldb] < zero) {
                        i32 nmi = n - i;
                        SLC_DSCAL(&nmi, &mone, &b[i + i * ldb], &ldb);
                    }
                }

                SLC_DCOPY(&nkh, &dwork[uiipt], &int1, &b[(kh + 1) * ldb + kl],
                         &ldb);
                if (kh > kl) {
                    SLC_DCOPY(&nkh, &dwork[uiipt + ldws], &int1,
                             &b[(kh + 1) * ldb + kh], &ldb);
                }
            }

            printf("\nDEBUG: Copying UI (%dx%d) to B[%d:%d,%d:%d]\n", kb, kb, kl, kl+kb-1, kl, kl+kb-1);
            printf("  UI matrix:\n");
            for (i32 ii = 0; ii < kb; ii++) {
                for (i32 jj = 0; jj < kb; jj++) {
                    printf("    UI[%d,%d] = %.16E\n", ii, jj, ui[ii + jj*2]);
                }
            }
            SLC_DLACPY("U", &kb, &kb, ui, &int2, &b[kl * ldb + kl], &ldb);

            printf("\nDEBUG: B diagonal after KB block processing and DLACPY:\n");
            for (i32 ii = 0; ii < n; ii++) {
                printf("  B[%d,%d] = %.16E\n", ii, ii, b[ii + ii*ldb]);
            }
            kh++;
        }

    } else {
        kl = n;

        while (kl > 0) {
            kh = kl - 1;
            if (kh == 0) {
                kl = 0;
                kb = 1;
            } else {
                if (kh > 0 && a[kh + (kh - 1) * lda] == zero) {
                    kl = kh;
                    kb = 1;
                } else {
                    kl = kh - 1;
                    kb = 2;
                }
            }
            kl1 = kl - 1;

            if (kb == 1) {
                delta1 = -a[kl * lda + kl];
                if (delta1 <= zero) {
                    *info = 3;
                    return;
                }
                delta1 = sqrt(delta1) * sqrt(e[kl * lde + kl]);
                t = (b[kl * ldb + kl] * smlnum) / sqtwo;
                if (t > delta1) {
                    scale1 = delta1 / t;
                    *scale = scale1 * (*scale);
                    SLC_DLASCL("Upper", &int0, &int0, &one, &scale1, &n, &n,
                              b, &ldb, &info1);
                }
                ui[0] = b[kl * ldb + kl] / delta1 / sqtwo;
                m1[0] = a[kl * lda + kl] / e[kl * lde + kl];
                m2[0] = (delta1 / e[kl * lde + kl]) * sqtwo;

            } else {
                sg03bx("C", "T", &a[kl * lda + kl], lda, &e[kl * lde + kl],
                      lde, &b[kl * ldb + kl], ldb, ui, int2, &scale1,
                      m1, int2, m2, int2, info);

                if (*info != 0) {
                    return;
                }

                if (scale1 != one) {
                    *scale = scale1 * (*scale);
                    SLC_DLASCL("Upper", &int0, &int0, &one, &scale1, &n, &n,
                              b, &ldb, &info1);
                }
            }

            if (kl > 0) {
                SLC_DGEMM("N", "T", &kl, &kb, &kb, &mone, &b[kl * ldb],
                         &ldb, m2, &int2, &zero, &dwork[uiipt], &ldws);

                SLC_DGEMM("N", "N", &kl, &kb, &kb, &mone, &a[kl * lda],
                         &lda, ui, &int2, &one, &dwork[uiipt], &ldws);

                SLC_DGEMM("N", "T", &kb, &kb, &kb, &one, ui, &int2, m1, &int2,
                         &zero, tm, &int2);

                SLC_DGEMM("N", "N", &kl, &kb, &kb, &mone, &e[kl * lde],
                         &lde, tm, &int2, &one, &dwork[uiipt], &ldws);

                SLC_DLASET("A", &kb, &kb, &zero, &one, tm, &int2);

                sg03bw("T", kl, kb, a, lda, tm, int2, e, lde, m1, int2,
                      &dwork[uiipt], ldws, &scale1, info);

                if (scale1 != one) {
                    *scale = scale1 * (*scale);
                    SLC_DLASCL("Upper", &int0, &int0, &one, &scale1, &n, &n,
                              b, &ldb, &info1);
                    i32 int4 = 4;
                    SLC_DSCAL(&int4, &scale1, ui, &int1);
                }

                SLC_DLACPY("A", &kl, &kb, &dwork[uiipt], &ldws,
                          &dwork[wpt], &ldws);

                SLC_DTRMM("L", "U", "N", "N", &kl, &kb, &one, e, &lde,
                         &dwork[wpt], &ldws);

                SLC_DGEMM("N", "N", &kl, &kb, &kb, &one, &e[kl * lde],
                         &lde, ui, &int2, &one, &dwork[wpt], &ldws);

                SLC_DLACPY("A", &kl, &kb, &b[kl * ldb], &ldb, &dwork[ypt],
                          &ldws);

                SLC_DGEMM("N", "N", &kl, &kb, &kb, &mone, &dwork[wpt],
                         &ldws, m2, &int2, &one, &dwork[ypt], &ldws);

                printf("\nDEBUG C: Starting QR factorization N-KH=%d, KB=%d\n", kl, kb);
                printf("  B diagonal before QR:\n");
                for (i32 ii = 0; ii < kl; ii++) {
                    printf("    B(%d,%d)=%.16e\n", ii, ii, b[ii + ii*ldb]);
                }

                l = ypt;
                for (j = 0; j < kb; j++) {
                    printf("\n--- QR ITERATION J=%d (processing column) ---\n", j);
                    for (i = kl - 1; i >= 0; i--) {
                        x = b[i + i * ldb];
                        t = dwork[l + i];
                        printf("  I=%d: X=B(%d,%d)=%.16e, T=DWORK(%ld)=%.16e\n",
                               kl - i - 1, i, i, x, l + i, t);
                        SLC_DLARTG(&x, &t, &c, &s, &r);
                        printf("       DLARTG output: R,C,S=%.16e %.16e %.16e\n", r, c, s);
                        b[i + i * ldb] = r;
                        printf("       Set B(%d,%d)=%.16e\n", i, i, r);
                        if (i > 0) {
                            printf("       Applying DROT to %d elements\n", i);
                            SLC_DROT(&i, b, &int1, &dwork[l], &int1, &c, &s);
                        }
                    }
                    printf("  After J=%d iteration, B diagonal:\n", j);
                    for (i32 ii = 0; ii < kl; ii++) {
                        printf("    B(%d,%d)=%.16e\n", ii, ii, b[ii + ii*ldb]);
                    }
                    l += ldws;
                }

                for (i = 0; i < kl; i++) {
                    if (b[i + i * ldb] < zero) {
                        i32 ip1 = i + 1;
                        SLC_DSCAL(&ip1, &mone, &b[i * ldb], &int1);
                    }
                }

                SLC_DLACPY("A", &kl, &kb, &dwork[uiipt], &ldws, &b[kl * ldb],
                          &ldb);
            }

            SLC_DLACPY("U", &kb, &kb, ui, &int2, &b[kl * ldb + kl], &ldb);
        }
    }
}
