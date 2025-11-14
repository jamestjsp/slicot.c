#include "slicot.h"
#include <math.h>
#include <stdbool.h>

static bool lsame(char ca, char cb) {
    if (ca >= 'a' && ca <= 'z') ca -= 32;
    if (cb >= 'a' && cb <= 'z') cb -= 32;
    return ca == cb;
}

static int32_t imin(int32_t a, int32_t b) {
    return a < b ? a : b;
}

static int32_t imax(int32_t a, int32_t b) {
    return a > b ? a : b;
}

void mb01qd(char type, int32_t m, int32_t n, int32_t kl, int32_t ku,
            double cfrom, double cto, int32_t nbl, const int32_t* nrows,
            double* a, int32_t lda, int32_t* info) {

    const double zero = 0.0, one = 1.0;
    int32_t itype;
    bool done, noblc;
    double smlnum, bignum, cfromc, ctoc, cfrom1, cto1, mul;
    int32_t i, j, jini, jfin, ifin, k, k1, k2, k3, k4;

    *info = 0;

    if (lsame(type, 'G')) {
        itype = 0;
    } else if (lsame(type, 'L')) {
        itype = 1;
    } else if (lsame(type, 'U')) {
        itype = 2;
    } else if (lsame(type, 'H')) {
        itype = 3;
    } else if (lsame(type, 'B')) {
        itype = 4;
    } else if (lsame(type, 'Q')) {
        itype = 5;
    } else {
        itype = 6;
    }

    if (imin(m, n) == 0)
        return;

    smlnum = 2.2250738585072014e-308;
    bignum = one / smlnum;

    cfromc = cfrom;
    ctoc = cto;

    do {
        cfrom1 = cfromc * smlnum;
        cto1 = ctoc / bignum;

        if (fabs(cfrom1) > fabs(ctoc) && ctoc != zero) {
            mul = smlnum;
            done = false;
            cfromc = cfrom1;
        } else if (fabs(cto1) > fabs(cfromc)) {
            mul = bignum;
            done = false;
            ctoc = cto1;
        } else {
            mul = ctoc / cfromc;
            done = true;
        }

        noblc = (nbl == 0);

        if (itype == 0) {
            for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                    a[j * lda + i] *= mul;
                }
            }
        } else if (itype == 1) {
            if (noblc) {
                for (j = 0; j < n; j++) {
                    for (i = j; i < m; i++) {
                        a[j * lda + i] *= mul;
                    }
                }
            } else {
                jfin = 0;
                for (k = 0; k < nbl; k++) {
                    jini = jfin;
                    jfin = jfin + nrows[k];
                    for (j = jini; j < jfin; j++) {
                        for (i = jini; i < m; i++) {
                            a[j * lda + i] *= mul;
                        }
                    }
                }
            }
        } else if (itype == 2) {
            if (noblc) {
                for (j = 0; j < n; j++) {
                    for (i = 0; i < imin(j + 1, m); i++) {
                        a[j * lda + i] *= mul;
                    }
                }
            } else {
                jfin = 0;
                for (k = 0; k < nbl; k++) {
                    jini = jfin;
                    jfin = jfin + nrows[k];
                    if (k == nbl - 1) jfin = n;
                    for (j = jini; j < jfin; j++) {
                        for (i = 0; i < imin(jfin, m); i++) {
                            a[j * lda + i] *= mul;
                        }
                    }
                }
            }
        } else if (itype == 3) {
            if (noblc) {
                for (j = 0; j < n; j++) {
                    for (i = 0; i < imin(j + 2, m); i++) {
                        a[j * lda + i] *= mul;
                    }
                }
            } else {
                jfin = 0;
                for (k = 0; k < nbl; k++) {
                    jini = jfin;
                    jfin = jfin + nrows[k];

                    if (k == nbl - 1) {
                        jfin = n;
                        ifin = n;
                    } else {
                        ifin = jfin + nrows[k + 1];
                    }

                    for (j = jini; j < jfin; j++) {
                        for (i = 0; i < imin(ifin, m); i++) {
                            a[j * lda + i] *= mul;
                        }
                    }
                }
            }
        } else if (itype == 4) {
            k3 = kl + 1;
            k4 = n + 1;
            for (j = 0; j < n; j++) {
                for (i = 0; i < imin(k3, k4 - j - 1); i++) {
                    a[j * lda + i] *= mul;
                }
            }
        } else if (itype == 5) {
            k1 = ku + 2;
            k3 = ku + 1;
            for (j = 0; j < n; j++) {
                for (i = imax(k1 - j - 2, 0); i < k3; i++) {
                    a[j * lda + i] *= mul;
                }
            }
        } else if (itype == 6) {
            k1 = kl + ku + 2;
            k2 = kl + 1;
            k3 = 2 * kl + ku + 1;
            k4 = kl + ku + 1 + m;
            for (j = 0; j < n; j++) {
                for (i = imax(k1 - j - 2, k2 - 1); i < imin(k3, k4 - j - 1); i++) {
                    a[j * lda + i] *= mul;
                }
            }
        }
    } while (!done);
}
