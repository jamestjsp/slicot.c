/*
 * SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 1996-2025, The SLICOT Team (original Fortran77 code)
 * Copyright (c) 2025, slicot.c contributors (C11 translation)
 */

/*
 * TB01XD - Special similarity transformation of dual state-space system
 *
 * Purpose:
 *   To apply a special transformation to a system given as a triple (A,B,C):
 *     A <-- P * A' * P,  B <-- P * C',  C <-- B' * P,
 *   where P is a matrix with 1 on the secondary diagonal, and with 0
 *   in the other entries. Matrix A can be specified as a band matrix.
 *   Optionally, matrix D of the system can be transposed.
 */

#include "slicot.h"
#include "slicot_blas.h"

void tb01xd(
    const char* jobd,
    const i32 n,
    const i32 m,
    const i32 p,
    const i32 kl,
    const i32 ku,
    f64* a,
    const i32 lda,
    f64* b,
    const i32 ldb,
    f64* c,
    const i32 ldc,
    f64* d,
    const i32 ldd,
    i32* info
)
{
    bool ljobd;
    i32 j, j1, lda1, maxmp, minmp, nm1;
    const i32 int1 = 1;
    i32 i, k;
    f64 temp;

    *info = 0;
    ljobd = (jobd[0] == 'D' || jobd[0] == 'd');
    maxmp = (m > p) ? m : p;
    minmp = (m < p) ? m : p;
    nm1 = n - 1;

    if (!ljobd && !(jobd[0] == 'Z' || jobd[0] == 'z')) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (p < 0) {
        *info = -4;
    } else if (kl < 0 || (nm1 >= 0 && kl > nm1) || (nm1 < 0 && kl > 0)) {
        *info = -5;
    } else if (ku < 0 || (nm1 >= 0 && ku > nm1) || (nm1 < 0 && ku > 0)) {
        *info = -6;
    } else if (lda < (n > 1 ? n : 1)) {
        *info = -8;
    } else if ((maxmp > 0 && ldb < (n > 1 ? n : 1)) || (minmp == 0 && ldb < 1)) {
        *info = -10;
    } else if (ldc < 1 || (n > 0 && ldc < maxmp)) {
        *info = -12;
    } else if (ldd < 1 || (ljobd && ldd < maxmp)) {
        *info = -14;
    }

    if (*info != 0) {
        return;
    }

    if (ljobd) {
        // Transpose D in-place
        // Fortran: DO J = 1, MAXMP
        for (j = 0; j < maxmp; j++) {
            if (j < minmp - 1) {
                // Swap below-diagonal with right-of-diagonal
                i32 count = minmp - j - 1;
                SLC_DSWAP(&count, &d[(j + 1) + j * ldd], &int1, &d[j + (j + 1) * ldd], &ldd);
            } else if (j >= p) {
                // Copy column j to row j (for j > p when m > p)
                SLC_DCOPY(&p, &d[j * ldd], &int1, &d[j], &ldd);
            } else if (j >= m) {
                // Copy row j to column j (for j > m when p > m)
                SLC_DCOPY(&m, &d[j], &ldd, &d[j * ldd], &int1);
            }
        }
    }

    if (n == 0) {
        return;
    }

    // Replace matrix A by P*A'*P (pertransposition)
    // Result[i,j] = A[n-1-j, n-1-i] (0-based indexing)
    // For band matrices, only process elements within the band

    (void)lda1; // Unused in simplified implementation

    // Process all elements within the band defined by kl and ku
    // For each (i,j), swap with (n-1-j, n-1-i) if within band
    for (i = 0; i < n; i++) {
        // Column range for row i within band: max(0, i-kl) to min(n-1, i+ku)
        i32 j_start = (i > kl) ? (i - kl) : 0;
        i32 j_end = ((i + ku) < (n - 1)) ? (i + ku) : (n - 1);

        for (j = j_start; j <= j_end; j++) {
            i32 i2 = n - 1 - j;
            i32 j2 = n - 1 - i;
            i32 idx1 = i + j * lda;
            i32 idx2 = i2 + j2 * lda;
            if (idx1 < idx2) {
                temp = a[idx1];
                a[idx1] = a[idx2];
                a[idx2] = temp;
            }
        }
    }

    // Replace matrix B by P*C' and matrix C by B'*P
    // P*C' means: take row j of C (= column j of C'), reverse it, put in column j of output B
    // B'*P means: take column j of B (= row j of B'), reverse it, put in row j of output C
    for (j = 0; j < maxmp; j++) {
        if (j < minmp) {
            // Swap B(:,j) with reversed C(j,:)
            // B(i,j) <-> C(j, n-1-i) for i = 0..n-1
            for (i = 0; i < n; i++) {
                i32 idx_b = i + j * ldb;
                i32 idx_c = j + (n - 1 - i) * ldc;
                temp = b[idx_b];
                b[idx_b] = c[idx_c];
                c[idx_c] = temp;
            }
        } else if (j >= p) {
            // Copy B(:,j) to reversed C(j,:)
            // j >= p means we're processing extra columns of B (since m > p)
            // C(j,:) is the j-th row of C (output), reversed <- B(:,j)
            for (i = 0; i < n; i++) {
                c[j + (n - 1 - i) * ldc] = b[i + j * ldb];
            }
        } else {
            // Copy reversed C(j,:) to B(:,j)
            // j >= m means we're processing extra rows of C (since p > m)
            // B(:,j) (output) <- reversed C(j,:)
            for (i = 0; i < n; i++) {
                b[i + j * ldb] = c[j + (n - 1 - i) * ldc];
            }
        }
    }
}
