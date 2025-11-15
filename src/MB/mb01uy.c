#include "slicot.h"
#include "slicot_blas.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

void mb01uy(
    const char* side, const char* uplo, const char* trans,
    const i32 m, const i32 n,
    const f64 alpha,
    f64* t, const i32 ldt,
    const f64* a, const i32 lda,
    f64* dwork, const i32 ldwork,
    i32* info
)
{
    const f64 zero = 0.0;

    bool lside, luplo, ltran;
    i32 k, mn, wrkmin;

    *info = 0;
    lside = (*side == 'L' || *side == 'l');
    luplo = (*uplo == 'U' || *uplo == 'u');
    ltran = (*trans == 'T' || *trans == 't' || *trans == 'C' || *trans == 'c');

    if (lside) {
        k = m;
    } else {
        k = n;
    }
    mn = (m < n) ? m : n;

    wrkmin = 1;
    if (alpha != zero && mn > 0) {
        wrkmin = (wrkmin > k) ? wrkmin : k;
    }

    if (ldwork == -1) {
        dwork[0] = (f64)(m * n);
        return;
    }

    if ((!lside && *side != 'R' && *side != 'r')) {
        *info = -1;
        return;
    }
    if ((!luplo && *uplo != 'L' && *uplo != 'l')) {
        *info = -2;
        return;
    }
    if ((!ltran && *trans != 'N' && *trans != 'n')) {
        *info = -3;
        return;
    }
    if (m < 0) {
        *info = -4;
        return;
    }
    if (n < 0) {
        *info = -5;
        return;
    }
    if (ldt < ((m > 1) ? m : 1) || (!lside && ldt < n)) {
        *info = -8;
        return;
    }
    if (lda < ((m > 1) ? m : 1)) {
        *info = -10;
        return;
    }
    if (ldwork < wrkmin) {
        dwork[0] = (f64)wrkmin;
        *info = -12;
        return;
    }

    if (mn == 0) {
        return;
    }

    if (alpha == zero) {
        SLC_DLASET("F", &m, &n, &zero, &zero, t, &ldt);
        return;
    }

    SLC_DLACPY("A", &m, &n, a, &lda, dwork, &m);

    SLC_DTRMM(side, uplo, trans, "N", &m, &n, &alpha, t, &ldt, dwork, &m);

    SLC_DLACPY("A", &m, &n, dwork, &m, t, &ldt);

    dwork[0] = (f64)(m * n);
}
