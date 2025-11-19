/*
 * SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 1996-2025, The SLICOT Team (original Fortran77 code)
 * Copyright (c) 2025, slicot.c contributors (C11 translation)
 */

#include "slicot.h"
#include "slicot_blas.h"
#include <math.h>

/**
 * @brief Calculate the output of a set of neural networks.
 *
 * Calculates the output of a set of neural networks with the structure:
 *
 *          - tanh(w1'*z+b1) -
 *        /      :             \
 *      z ---    :           --- sum(ws(i)*...)+ b(n+1)  --- y,
 *        \      :             /
 *          - tanh(wn'*z+bn) -
 *
 * given the input z and the parameter vectors wi, ws, and b.
 *
 * @param[in] nsmp The number of training samples. NSMP >= 0.
 * @param[in] nz The length of each input sample. NZ >= 0.
 * @param[in] l The length of each output sample. L >= 0.
 * @param[in] ipar Integer parameters. ipar[0] must contain NN (neurons).
 * @param[in] lipar Length of ipar.
 * @param[in] wb Weights and biases vector.
 * @param[in] lwb Length of wb.
 * @param[in] z Input samples matrix (nsmp x nz).
 * @param[in] ldz Leading dimension of z.
 * @param[out] y Output samples matrix (nsmp x l).
 * @param[in] ldy Leading dimension of y.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01ay(i32 nsmp, i32 nz, i32 l, const i32 *ipar, i32 lipar, 
            const f64 *wb, i32 lwb, const f64 *z, i32 ldz, 
            f64 *y, i32 ldy, f64 *dwork, i32 ldwork, i32 *info)
{
    /* Local variables */
    i32 i, ib, j, k, ldwb, lj, lk, m, mf, nn, nv, ws;
    f64 bignum, df, smlnum, tmp;
    bool last;
    
    /* Constants */
    f64 zero = 0.0;
    f64 one = 1.0;
    f64 two = 2.0;
    f64 neg_two = -2.0;

    *info = 0;
    nn = ipar[0];
    ldwb = nn * (nz + 2) + 1;

    /* Argument checks */
    if (nsmp < 0) *info = -1;
    else if (nz < 0) *info = -2;
    else if (l < 0) *info = -3;
    else if (nn < 0) *info = -4;
    else if (lipar < 1) *info = -5;
    else if (lwb < ldwb * l) *info = -7;
    else if (ldz < 1 || (nsmp > 0 && ldz < nsmp)) *info = -9;
    else if (ldy < 1 || (nsmp > 0 && ldy < nsmp)) *info = -11;
    else if (ldwork < 2 * nn) *info = -13;

    if (*info != 0) {
        i32 err_code = -(*info);
        SLC_XERBLA("NF01AY", &err_code);
        return;
    }

    if (nsmp == 0 || l == 0) return;

    /* Set parameters to avoid overflows */
    smlnum = SLC_DLAMCH("Safe minimum") / SLC_DLAMCH("Precision");
    bignum = one / smlnum;
    SLC_DLABAD(&smlnum, &bignum);
    smlnum = log(smlnum);
    bignum = log(bignum);

    ws = nz * nn; /* 0-based offset for ws part (starts after nz*nn elements) */
    ib = ws + nn - 1; /* This seems wrong compared to Fortran */
    
    /* Fortran: WS = NZ*NN + 1
       IB = WS + NN - 1
       WB(IB+NN+1+LK)
       
       Let's map to 0-based offsets.
       W part: 0 to NZ*NN - 1.
       WS part: starts at NZ*NN. Length NN.
       So WS_offset = NZ*NN.
       B part: starts at NZ*NN + NN. Length NN+1.
       
       Fortran WS = NZ*NN + 1.
       Fortran IB = WS + NN - 1 = NZ*NN + NN.
       So IB is the index of the last element of WS?
       
       Loop J=1,NN: DWORK(J) = TWO * WB(IB+J+LK)
       If J=1, WB(IB+1+LK).
       If IB = NZ*NN + NN. IB+1 = NZ*NN + NN + 1.
       This is the first element of B. Correct.
       
       So in C:
       ws_offset = nz * nn;
       b_offset = ws_offset + nn;
       
       WB[b_offset + j + lk] corresponds to b(j+1).
    */
    
    ws = nz * nn;
    ib = ws + nn; /* Points to start of B vector in 0-based indexing */
    
    lk = 0;
    
    if (nz == 0 || nn == 0) {
        nv = 2;
    } else {
        nv = (ldwork - nn) / nn;
    }

    if (nv > 2) {
        mf = (nsmp / nv) * nv;
        last = (nsmp % nv) != 0;

        /* BLAS 3 calculations */
        for (k = 0; k < l; k++) {
            /* b(n+1) is at b_offset + nn */
            tmp = wb[ib + nn + lk];

            for (j = 0; j < nn; j++) {
                dwork[j] = two * wb[ib + j + lk];
            }

            for (i = 0; i < mf; i += nv) {
                /* Compute -2*[w1 ... wn]' * Z' */
                /* WB(1+LK) -> &wb[lk]
                   Z(I,1) -> &z[i] (start of block of rows)
                   DWORK(NN+1) -> &dwork[nn]
                */
                SLC_DGEMM("Transpose", "Transpose", &nn, &nv, &nz, &neg_two,
                          &wb[lk], &nz, &z[i], &ldz, &zero, &dwork[nn], &nn);

                lj = nn;
                for (m = 0; m < nv; m++) {
                    for (j = 0; j < nn; j++) {
                        lj++; /* dwork[lj] is result element */
                        /* Using lj-1 because of 0-based indexing vs Fortran loop */
                        /* Wait. In Fortran:
                           LJ = NN
                           DO M=1,NV
                             DO J=1,NN
                               LJ = LJ + 1
                               DF = DWORK(LJ) - DWORK(J)
                           
                           If M=1, J=1: LJ = NN+1.
                           DWORK(NN+1) is first element of result matrix C.
                           
                           In C: dwork starts at 0.
                           dwork[nn] is start of result C.
                           So dwork[nn + m*nn + j] is element (j, m) of C.
                           Since C is NN x NV, stored column major (leading dim NN).
                           So element (j, m) is indeed at index j + m*nn.
                           
                           Let's simplify loop:
                        */
                        i32 idx = nn + m * nn + j;
                        df = dwork[idx] - dwork[j];
                        
                        if (fabs(df) >= bignum) {
                            if (df > zero) dwork[idx] = -one;
                            else dwork[idx] = one;
                        } else if (fabs(df) <= smlnum) {
                            dwork[idx] = zero;
                        } else {
                            dwork[idx] = two / (one + exp(df)) - one;
                        }
                    }
                }

                /* Y(I, K+1) = TMP */
                /* Y is column major. Element (I, K) is Y[I + K*LDY].
                   Wait, Fortran Y(I, K+1). 1-based.
                   0-based: Y[i + k*ldy].
                */
                y[i + k * ldy] = tmp;
                
                /* CALL DCOPY( NV-1, Y(I, K+1), 0, Y(I+1, K+1), 1 ) */
                /* Copies TMP to next NV-1 elements of Y column */
                i32 nv_minus_1 = nv - 1;
                i32 inc_0 = 0;
                i32 inc_1 = 1;
                if (nv_minus_1 > 0) {
                    SLC_DCOPY(&nv_minus_1, &y[i + k * ldy], &inc_0, 
                              &y[(i + 1) + k * ldy], &inc_1);
                }

                /* CALL DGEMV( 'Transpose', NN, NV, ONE, DWORK(NN+1), NN,
                               WB(WS+LK), 1, ONE, Y(I, K+1), 1 ) */
                /* WB(WS+LK) -> &wb[ws + lk] (weights ws)
                   DWORK(NN+1) -> &dwork[nn] (tanh outputs)
                   
                   Y(I, K+1) -> &y[i + k*ldy]. Stride 1.
                   
                   DWORK(NN+1) is NN x NV matrix of tanh values.
                   We want: y_vec = y_vec + DWORK' * ws
                   
                   DGEMV: y = alpha*A'*x + beta*y
                   A = DWORK(NN+1). NN x NV.
                   x = ws. Length NN.
                   y = subvector of Y. Length NV.
                   
                   So Y[i...i+nv-1, k] += (NN x NV)' * (NN x 1) = (NV x 1).
                   Correct.
                */
                SLC_DGEMV("Transpose", &nn, &nv, &one, &dwork[nn], &nn,
                          &wb[ws + lk], &inc_1, &one, &y[i + k * ldy], &inc_1);
            }

            if (last) {
                nv = nsmp - mf;
                i = mf; /* Start index */
                
                SLC_DGEMM("Transpose", "Transpose", &nn, &nv, &nz, &neg_two,
                          &wb[lk], &nz, &z[i], &ldz, &zero, &dwork[nn], &nn);
                
                for (m = 0; m < nv; m++) {
                    for (j = 0; j < nn; j++) {
                        i32 idx = nn + m * nn + j;
                        df = dwork[idx] - dwork[j];
                        if (fabs(df) >= bignum) {
                            if (df > zero) dwork[idx] = -one;
                            else dwork[idx] = one;
                        } else if (fabs(df) <= smlnum) {
                            dwork[idx] = zero;
                        } else {
                            dwork[idx] = two / (one + exp(df)) - one;
                        }
                    }
                }
                
                y[i + k * ldy] = tmp;
                i32 nv_minus_1 = nv - 1;
                i32 inc_0 = 0;
                i32 inc_1 = 1;
                if (nv_minus_1 > 0) {
                    SLC_DCOPY(&nv_minus_1, &y[i + k * ldy], &inc_0, 
                              &y[(i + 1) + k * ldy], &inc_1);
                }
                
                SLC_DGEMV("Transpose", &nn, &nv, &one, &dwork[nn], &nn,
                          &wb[ws + lk], &inc_1, &one, &y[i + k * ldy], &inc_1);
            }
            
            lk += ldwb;
        }

    } else {
        /* BLAS 2 calculations only */
        for (k = 0; k < l; k++) {
            tmp = wb[ib + nn + lk];

            for (j = 0; j < nn; j++) {
                dwork[j] = two * wb[ib + j + lk];
            }

            for (i = 0; i < nsmp; i++) {
                if (nz == 0) {
                    for (j = 0; j < nn; j++) dwork[nn + j] = zero;
                } else {
                    i32 inc_1 = 1;
                    SLC_DGEMV("Transpose", &nz, &nn, &neg_two, &wb[lk], &nz,
                              &z[i], &ldz, &zero, &dwork[nn], &inc_1);
                }

                for (j = nn; j < 2 * nn; j++) {
                    /* j goes from nn to 2*nn-1.
                       DWORK(J) in Fortran (1-based) maps to dwork[j-1] in C?
                       Wait, DWORK indices in C are 0-based.
                       Fortran: DO 90 J = NN + 1, 2*NN.
                       DWORK(J).
                       Corresponding C index: j.
                       
                       DF = DWORK(J) - DWORK(J-NN).
                       Note: DWORK[0..NN-1] holds bias values.
                       DWORK[NN..2*NN-1] holds matrix product results.
                       
                       So indices match if we use 0-based everywhere.
                       Loop j from nn to 2*nn - 1.
                       df = dwork[j] - dwork[j - nn].
                    */
                    df = dwork[j] - dwork[j - nn];
                    
                    if (fabs(df) >= bignum) {
                        if (df > zero) dwork[j] = -one;
                        else dwork[j] = one;
                    } else if (fabs(df) <= smlnum) {
                        dwork[j] = zero;
                    } else {
                        dwork[j] = two / (one + exp(df)) - one;
                    }
                }
                
                i32 inc_1 = 1;
                y[i + k * ldy] = SLC_DDOT(&nn, &wb[ws + lk], &inc_1, &dwork[nn], &inc_1) + tmp;
            }
            lk += ldwb;
        }
    }
}
