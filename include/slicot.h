/*
 * SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 1996-2025, The SLICOT Team (original Fortran77 code)
 * Copyright (c) 2025, slicot.c contributors (C11 translation)
 */

#ifndef SLICOT_H
#define SLICOT_H

#include "slicot_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Main header for SLICOT C library */

/**
 * @brief Multiply matrix by scalar CTO/CFROM without overflow/underflow.
 *
 * Multiplies M-by-N real matrix A by scalar CTO/CFROM, preventing overflow and
 * underflow through iterative scaling. The algorithm decomposes the scaling
 * factor into products of safe intermediate values (DBL_MIN or DBL_MAX), applying
 * them sequentially until the final result is achieved. This ensures that
 * intermediate values remain representable even when the final result is near
 * machine limits.
 *
 * Supports full, triangular, Hessenberg, and banded matrix storage with optional
 * block structure for efficient handling of structured matrices.
 *
 * Based on LAPACK routine DLASCL with extensions for block-structured matrices.
 * For efficiency, input parameters are not validated (caller responsibility).
 *
 * @param[in] type Matrix storage type:
 *                 'G' = full matrix
 *                 'L' = (block) lower triangular
 *                 'U' = (block) upper triangular
 *                 'H' = (block) upper Hessenberg
 *                 'B' = symmetric band (lower half)
 *                 'Q' = symmetric band (upper half)
 *                 'Z' = general band
 * @param[in] m Number of rows (M >= 0)
 * @param[in] n Number of columns (N >= 0)
 * @param[in] kl Lower bandwidth (for 'B', 'Q', 'Z' types)
 * @param[in] ku Upper bandwidth (for 'B', 'Q', 'Z' types)
 * @param[in] cfrom Denominator scalar (caller must ensure nonzero)
 * @param[in] cto Numerator scalar
 * @param[in] nbl Number of diagonal blocks (0 = no block structure)
 * @param[in] nrows Block sizes, dimension max(1,nbl), may be NULL if nbl=0
 * @param[in,out] a Matrix array, dimension (lda,n), column-major storage
 * @param[in] lda Leading dimension (lda >= max(1,m))
 * @param[out] info Exit code (always 0, reserved for future use)
 */
void mb01qd(char type, i32 m, i32 n, i32 kl, i32 ku,
            f64 cfrom, f64 cto, i32 nbl, const i32* nrows,
            f64* a, i32 lda, i32* info);

/**
 * @brief Bilinear transformation of state-space system.
 *
 * Performs discrete-time <-> continuous-time conversion via bilinear
 * transformation of the state-space matrices (A,B,C,D).
 *
 * For TYPE='D' (discrete -> continuous):
 *   A_out = beta * (alpha*I + A)^{-1} * (A - alpha*I)
 *   B_out = sqrt(2*alpha*beta) * (alpha*I + A)^{-1} * B
 *   C_out = sqrt(2*alpha*beta) * C * (alpha*I + A)^{-1}
 *   D_out = D - C * (alpha*I + A)^{-1} * B
 *
 * For TYPE='C' (continuous -> discrete):
 *   A_out = alpha * (beta*I - A)^{-1} * (beta*I + A)
 *   B_out = sqrt(2*alpha*beta) * (beta*I - A)^{-1} * B
 *   C_out = sqrt(2*alpha*beta) * C * (beta*I - A)^{-1}
 *   D_out = D + C * (beta*I - A)^{-1} * B
 *
 * @param[in] type 'D' for discrete->continuous, 'C' for continuous->discrete
 * @param[in] n Order of state matrix A (n >= 0)
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] p Number of system outputs (p >= 0)
 * @param[in] alpha Bilinear transformation parameter (alpha != 0)
 * @param[in] beta Bilinear transformation parameter (beta != 0)
 * @param[in,out] a State matrix, dimension (lda, n). On exit, transformed A.
 * @param[in] lda Leading dimension of a (lda >= max(1, n))
 * @param[in,out] b Input matrix, dimension (ldb, m). On exit, transformed B.
 * @param[in] ldb Leading dimension of b (ldb >= max(1, n))
 * @param[in,out] c Output matrix, dimension (ldc, n). On exit, transformed C.
 * @param[in] ldc Leading dimension of c (ldc >= max(1, p))
 * @param[in,out] d Feedthrough matrix, dimension (ldd, m). On exit, transformed D.
 * @param[in] ldd Leading dimension of d (ldd >= max(1, p))
 * @param[out] iwork Integer workspace, dimension (n)
 * @param[out] dwork Double workspace, dimension (ldwork)
 * @param[in] ldwork Workspace size (ldwork >= max(1, n))
 * @return 0 on success, -i if parameter i is invalid, 1 if (alpha*I + A) singular,
 *         2 if (beta*I - A) singular
 */
i32 ab04md(char type, i32 n, i32 m, i32 p, f64 alpha, f64 beta,
           f64* a, i32 lda, f64* b, i32 ldb, f64* c, i32 ldc, f64* d, i32 ldd,
           i32* iwork, f64* dwork, i32 ldwork);

/**
 * @brief Cascade (series) inter-connection of two state-space systems.
 *
 * Computes the state-space model (A,B,C,D) for the cascaded connection
 * of two systems G1 and G2, where output of G1 feeds input of G2:
 * Y = G2(G1(U))
 *
 * For UPLO='L' (lower block diagonal):
 *   A = [A1,    0   ]    B = [ B1   ]    C = [D2*C1, C2]    D = D2*D1
 *       [B2*C1, A2  ]        [B2*D1 ]
 *
 * For UPLO='U' (upper block diagonal):
 *   A = [A2,    B2*C1]    B = [B2*D1]    C = [C2, D2*C1]    D = D2*D1
 *       [0,     A1   ]        [ B1  ]
 *
 * @param[in] uplo 'L' for lower block diagonal, 'U' for upper block diagonal
 * @param[in] over 'N' no overlap, 'O' overlap arrays (workspace required)
 * @param[in] n1 Number of states in first system (n1 >= 0)
 * @param[in] m1 Number of inputs to first system (m1 >= 0)
 * @param[in] p1 Number of outputs from G1 = inputs to G2 (p1 >= 0)
 * @param[in] n2 Number of states in second system (n2 >= 0)
 * @param[in] p2 Number of outputs from second system (p2 >= 0)
 * @param[in] a1 State matrix of G1, dimension (lda1, n1)
 * @param[in] lda1 Leading dimension of a1 (lda1 >= max(1, n1))
 * @param[in] b1 Input matrix of G1, dimension (ldb1, m1)
 * @param[in] ldb1 Leading dimension of b1 (ldb1 >= max(1, n1))
 * @param[in] c1 Output matrix of G1, dimension (ldc1, n1)
 * @param[in] ldc1 Leading dimension of c1 (ldc1 >= max(1, p1) if n1 > 0)
 * @param[in] d1 Feedthrough matrix of G1, dimension (ldd1, m1)
 * @param[in] ldd1 Leading dimension of d1 (ldd1 >= max(1, p1))
 * @param[in] a2 State matrix of G2, dimension (lda2, n2)
 * @param[in] lda2 Leading dimension of a2 (lda2 >= max(1, n2))
 * @param[in] b2 Input matrix of G2, dimension (ldb2, p1)
 * @param[in] ldb2 Leading dimension of b2 (ldb2 >= max(1, n2))
 * @param[in] c2 Output matrix of G2, dimension (ldc2, n2)
 * @param[in] ldc2 Leading dimension of c2 (ldc2 >= max(1, p2) if n2 > 0)
 * @param[in] d2 Feedthrough matrix of G2, dimension (ldd2, p1)
 * @param[in] ldd2 Leading dimension of d2 (ldd2 >= max(1, p2))
 * @param[out] n Total state order of cascaded system (n = n1 + n2)
 * @param[out] a State matrix of cascaded system, dimension (lda, n1+n2)
 * @param[in] lda Leading dimension of a (lda >= max(1, n1+n2))
 * @param[out] b Input matrix of cascaded system, dimension (ldb, m1)
 * @param[in] ldb Leading dimension of b (ldb >= max(1, n1+n2))
 * @param[out] c Output matrix of cascaded system, dimension (ldc, n1+n2)
 * @param[in] ldc Leading dimension of c (ldc >= max(1, p2) if n1+n2 > 0)
 * @param[out] d Feedthrough matrix of cascaded system, dimension (ldd, m1)
 * @param[in] ldd Leading dimension of d (ldd >= max(1, p2))
 * @param[out] dwork Workspace, dimension (ldwork)
 * @param[in] ldwork Workspace size (ldwork >= max(1, p1*max(n1,m1,n2,p2)) if over='O')
 * @return 0 on success, -i if parameter i is invalid
 */
i32 ab05md(char uplo, char over, i32 n1, i32 m1, i32 p1, i32 n2, i32 p2,
           const f64* a1, i32 lda1, const f64* b1, i32 ldb1,
           const f64* c1, i32 ldc1, const f64* d1, i32 ldd1,
           const f64* a2, i32 lda2, const f64* b2, i32 ldb2,
           const f64* c2, i32 ldc2, const f64* d2, i32 ldd2,
           i32* n, f64* a, i32 lda, f64* b, i32 ldb,
           f64* c, i32 ldc, f64* d, i32 ldd,
           f64* dwork, i32 ldwork);

/**
 * @brief Feedback inter-connection of two state-space systems.
 *
 * Computes the state-space model (A,B,C,D) for the feedback connection
 * of two systems G1 and G2:
 *   U = U1 + alpha*Y2,  Y = Y1 = U2
 *   alpha = +1: positive feedback
 *   alpha = -1: negative feedback
 *
 * The interconnection matrices are:
 *   E21 = (I + alpha*D1*D2)^-1
 *   E12 = I - alpha*D2*E21*D1
 *
 * Matrix A:
 *   [A1 - alpha*B1*E12*D2*C1,    -alpha*B1*E12*C2    ]
 *   [B2*E21*C1,                   A2 - alpha*B2*E21*D1*C2]
 *
 * Matrix B:  [B1*E12; B2*E21*D1]
 * Matrix C:  [E21*C1, -alpha*E21*D1*C2]
 * Matrix D:  E21*D1
 *
 * @param[in] over 'N' no overlap, 'O' overlap arrays (workspace required)
 * @param[in] n1 Number of states in first system (n1 >= 0)
 * @param[in] m1 Number of inputs to first system and outputs from G2 (m1 >= 0)
 * @param[in] p1 Number of outputs from G1 and inputs to G2 (p1 >= 0)
 * @param[in] n2 Number of states in second system (n2 >= 0)
 * @param[in] alpha Feedback coefficient (+1 positive, -1 negative feedback)
 * @param[in] a1 State matrix of G1, dimension (lda1, n1)
 * @param[in] lda1 Leading dimension of a1 (lda1 >= max(1, n1))
 * @param[in] b1 Input matrix of G1, dimension (ldb1, m1)
 * @param[in] ldb1 Leading dimension of b1 (ldb1 >= max(1, n1))
 * @param[in] c1 Output matrix of G1, dimension (ldc1, n1)
 * @param[in] ldc1 Leading dimension of c1 (ldc1 >= max(1, p1) if n1 > 0)
 * @param[in] d1 Feedthrough matrix of G1, dimension (ldd1, m1)
 * @param[in] ldd1 Leading dimension of d1 (ldd1 >= max(1, p1))
 * @param[in] a2 State matrix of G2, dimension (lda2, n2)
 * @param[in] lda2 Leading dimension of a2 (lda2 >= max(1, n2))
 * @param[in] b2 Input matrix of G2, dimension (ldb2, p1)
 * @param[in] ldb2 Leading dimension of b2 (ldb2 >= max(1, n2))
 * @param[in] c2 Output matrix of G2, dimension (ldc2, n2)
 * @param[in] ldc2 Leading dimension of c2 (ldc2 >= max(1, m1) if n2 > 0)
 * @param[in] d2 Feedthrough matrix of G2, dimension (ldd2, p1)
 * @param[in] ldd2 Leading dimension of d2 (ldd2 >= max(1, m1))
 * @param[out] n Total state order (n = n1 + n2)
 * @param[out] a State matrix, dimension (lda, n1+n2)
 * @param[in] lda Leading dimension of a (lda >= max(1, n1+n2))
 * @param[out] b Input matrix, dimension (ldb, m1)
 * @param[in] ldb Leading dimension of b (ldb >= max(1, n1+n2))
 * @param[out] c Output matrix, dimension (ldc, n1+n2)
 * @param[in] ldc Leading dimension of c (ldc >= max(1, p1) if n1+n2 > 0)
 * @param[out] d Feedthrough matrix, dimension (ldd, m1)
 * @param[in] ldd Leading dimension of d (ldd >= max(1, p1))
 * @param[out] iwork Integer workspace, dimension (p1)
 * @param[out] dwork Double workspace, dimension (ldwork)
 * @param[in] ldwork Workspace size
 * @return 0 on success, -i if parameter i is invalid, >0 if singular
 */
i32 ab05nd(char over, i32 n1, i32 m1, i32 p1, i32 n2, f64 alpha,
           const f64* a1, i32 lda1, const f64* b1, i32 ldb1,
           const f64* c1, i32 ldc1, const f64* d1, i32 ldd1,
           const f64* a2, i32 lda2, const f64* b2, i32 ldb2,
           const f64* c2, i32 ldc2, const f64* d2, i32 ldd2,
           i32* n, f64* a, i32 lda, f64* b, i32 ldb,
           f64* c, i32 ldc, f64* d, i32 ldd,
           i32* iwork, f64* dwork, i32 ldwork);

/**
 * @brief Rowwise concatenation of two state-space systems.
 *
 * Computes the state-space model (A,B,C,D) for rowwise concatenation
 * (parallel inter-connection on outputs, with separate inputs) of two
 * systems G1 and G2:
 *   Y = G1*U1 + alpha*G2*U2
 *
 * The combined system has:
 *   A = [[A1, 0], [0, A2]]     (block diagonal)
 *   B = [[B1, 0], [0, B2]]     (block diagonal)
 *   C = [C1, alpha*C2]         (rowwise concatenation)
 *   D = [D1, alpha*D2]         (rowwise concatenation)
 *
 * @param[in] over 'N' no overlap, 'O' overlap arrays A1/A, B1/B, C1/C, D1/D
 * @param[in] n1 Number of states in first system (n1 >= 0)
 * @param[in] m1 Number of inputs to first system (m1 >= 0)
 * @param[in] p1 Number of outputs from each system (p1 >= 0)
 * @param[in] n2 Number of states in second system (n2 >= 0)
 * @param[in] m2 Number of inputs to second system (m2 >= 0)
 * @param[in] alpha Coefficient multiplying second system output
 * @param[in] a1 State matrix of G1, dimension (lda1, n1)
 * @param[in] lda1 Leading dimension of a1 (lda1 >= max(1, n1))
 * @param[in] b1 Input matrix of G1, dimension (ldb1, m1)
 * @param[in] ldb1 Leading dimension of b1 (ldb1 >= max(1, n1))
 * @param[in] c1 Output matrix of G1, dimension (ldc1, n1)
 * @param[in] ldc1 Leading dimension of c1 (ldc1 >= max(1, p1) if n1 > 0)
 * @param[in] d1 Feedthrough matrix of G1, dimension (ldd1, m1)
 * @param[in] ldd1 Leading dimension of d1 (ldd1 >= max(1, p1))
 * @param[in] a2 State matrix of G2, dimension (lda2, n2)
 * @param[in] lda2 Leading dimension of a2 (lda2 >= max(1, n2))
 * @param[in] b2 Input matrix of G2, dimension (ldb2, m2)
 * @param[in] ldb2 Leading dimension of b2 (ldb2 >= max(1, n2))
 * @param[in] c2 Output matrix of G2, dimension (ldc2, n2)
 * @param[in] ldc2 Leading dimension of c2 (ldc2 >= max(1, p1) if n2 > 0)
 * @param[in] d2 Feedthrough matrix of G2, dimension (ldd2, m2)
 * @param[in] ldd2 Leading dimension of d2 (ldd2 >= max(1, p1))
 * @param[out] n Total state order (n = n1 + n2)
 * @param[out] m Total inputs (m = m1 + m2)
 * @param[out] a State matrix, dimension (lda, n1+n2)
 * @param[in] lda Leading dimension of a (lda >= max(1, n1+n2))
 * @param[out] b Input matrix, dimension (ldb, m1+m2)
 * @param[in] ldb Leading dimension of b (ldb >= max(1, n1+n2))
 * @param[out] c Output matrix, dimension (ldc, n1+n2)
 * @param[in] ldc Leading dimension of c (ldc >= max(1, p1) if n1+n2 > 0)
 * @param[out] d Feedthrough matrix, dimension (ldd, m1+m2)
 * @param[in] ldd Leading dimension of d (ldd >= max(1, p1))
 * @return 0 on success, -i if parameter i is invalid
 */
i32 ab05od(char over, i32 n1, i32 m1, i32 p1, i32 n2, i32 m2, f64 alpha,
           const f64* a1, i32 lda1, const f64* b1, i32 ldb1,
           const f64* c1, i32 ldc1, const f64* d1, i32 ldd1,
           const f64* a2, i32 lda2, const f64* b2, i32 ldb2,
           const f64* c2, i32 ldc2, const f64* d2, i32 ldd2,
           i32* n, i32* m, f64* a, i32 lda, f64* b, i32 ldb,
           f64* c, i32 ldc, f64* d, i32 ldd);

/**
 * @brief Compute the inverse of a linear system.
 *
 * Computes the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D):
 *   Ai = A - B*D^-1*C,  Bi = -B*D^-1,  Ci = D^-1*C,  Di = D^-1.
 *
 * @param[in] n Order of state matrix A (n >= 0)
 * @param[in] m Number of system inputs/outputs (m >= 0, square system)
 * @param[in,out] a State matrix, dimension (lda, n). On exit, Ai.
 * @param[in] lda Leading dimension of a (lda >= max(1, n))
 * @param[in,out] b Input matrix, dimension (ldb, m). On exit, Bi.
 * @param[in] ldb Leading dimension of b (ldb >= max(1, n))
 * @param[in,out] c Output matrix, dimension (ldc, n). On exit, Ci.
 * @param[in] ldc Leading dimension of c (ldc >= max(1, m))
 * @param[in,out] d Feedthrough matrix, dimension (ldd, m). On exit, Di.
 * @param[in] ldd Leading dimension of d (ldd >= max(1, m))
 * @param[out] rcond Reciprocal condition number of D
 * @param[out] iwork Integer workspace, dimension (2*m)
 * @param[out] dwork Double workspace, dimension (ldwork)
 * @param[in] ldwork Workspace size (ldwork >= max(1, 4*m))
 * @return 0 on success, -i if parameter i is invalid, 1..m if D is singular
 *         (i-th diagonal element zero), m+1 if D is numerically singular
 */
i32 ab07nd(i32 n, i32 m, f64* a, i32 lda, f64* b, i32 ldb,
           f64* c, i32 ldc, f64* d, i32 ldd, f64* rcond,
           i32* iwork, f64* dwork, i32 ldwork);

/**
 * @brief Compute matrix product T := alpha*op(T)*A or T := alpha*A*op(T).
 *
 * Computes one of the matrix products:
 *   T := alpha*op(T)*A (SIDE='L'), or
 *   T := alpha*A*op(T) (SIDE='R'),
 * where alpha is a scalar, A is M-by-N, T is triangular, and op(T) is T or T^T.
 * Uses block algorithm with BLAS 3 when possible. Result overwrites T.
 *
 * @param[in] side 'L' for T:=alpha*op(T)*A, 'R' for T:=alpha*A*op(T)
 * @param[in] uplo 'U' for upper triangular T, 'L' for lower triangular
 * @param[in] trans 'N' for op(T)=T, 'T'/'C' for op(T)=T^T
 * @param[in] m Number of rows of A (m >= 0)
 * @param[in] n Number of columns of A (n >= 0)
 * @param[in] alpha Scalar multiplier (alpha=0: T set to zero)
 * @param[in,out] t DOUBLE PRECISION array, dimension (ldt,max(K,N)) if SIDE='L',
 *                  (ldt,K) if SIDE='R', where K=M if SIDE='L', K=N if SIDE='R'.
 *                  In: K-by-K triangular matrix T
 *                  Out: M-by-N result matrix
 * @param[in] ldt Leading dimension of T (ldt >= max(1,M) if SIDE='L',
 *                ldt >= max(1,M,N) if SIDE='R')
 * @param[in] a DOUBLE PRECISION array, dimension (lda,n), M-by-N matrix A
 * @param[in] lda Leading dimension of A (lda >= max(1,M))
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork), workspace
 *                   On exit, dwork[0] returns optimal ldwork
 * @param[in] ldwork Workspace size (ldwork >= 1 if alpha=0 or min(M,N)=0,
 *                   ldwork >= M if SIDE='L', ldwork >= N if SIDE='R'.
 *                   If ldwork=-1, workspace query mode)
 * @param[out] info Exit code (0=success, <0=invalid parameter)
 */
void mb01uy(
    const char* side, const char* uplo, const char* trans,
    const i32 m, const i32 n,
    const f64 alpha,
    f64* t, const i32 ldt,
    const f64* a, const i32 lda,
    f64* dwork, const i32 ldwork,
    i32* info
);

/**
 * @brief Incremental rank estimation for QR factorization.
 *
 * Computes (optionally) a rank-revealing QR factorization of real M-by-N matrix A
 * and estimates effective rank using incremental condition estimation.
 *
 * Uses QR with column pivoting: A*P = Q*R, where R = [R11 R12; 0 R22].
 * R11 is largest leading submatrix with estimated condition < 1/RCOND.
 * Order of R11 (RANK) is effective rank of A.
 *
 * @param[in] jobqr 'Q' = perform QR factorization, 'N' = use existing R in A
 * @param[in] m Number of rows of A (m >= 0)
 * @param[in] n Number of columns of A (n >= 0)
 * @param[in,out] a DOUBLE PRECISION array, dimension (lda,n)
 *                  If JOBQR='Q': In: M-by-N matrix A, Out: QR factorization
 *                  If JOBQR='N': In/Out: Upper triangular R from QR
 * @param[in] lda Leading dimension of A (lda >= max(1,m))
 * @param[in,out] jpvt INTEGER array, dimension (n)
 *                     If JOBQR='Q': In: initial column flags (0=free),
 *                                   Out: pivot permutation
 *                     If JOBQR='N': not referenced
 * @param[in] rcond Rank threshold: condition < 1/RCOND (rcond >= 0)
 * @param[in] svlmax Largest singular value estimate of parent matrix (svlmax >= 0)
 *                   Use 0 if A is standalone
 * @param[out] tau DOUBLE PRECISION array, dimension (min(m,n))
 *                 Scalar factors of elementary reflectors (if JOBQR='Q')
 * @param[out] rank Effective rank of A
 * @param[out] sval DOUBLE PRECISION array, dimension (3)
 *                  sval[0]: largest singular value of R(1:rank,1:rank)
 *                  sval[1]: smallest singular value of R(1:rank,1:rank)
 *                  sval[2]: smallest singular value of R(1:rank+1,1:rank+1)
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   On exit, dwork[0] = optimal ldwork
 * @param[in] ldwork Length of dwork
 *                   JOBQR='Q': ldwork >= 3*n+1 (prefer 2*n+(n+1)*NB)
 *                   JOBQR='N': ldwork >= max(1,2*min(m,n))
 *                   If ldwork=-1, workspace query
 * @param[out] info Exit code (0=success, <0=invalid parameter -info)
 */
void mb03od(
    const char* jobqr,
    const i32 m,
    const i32 n,
    f64* a,
    const i32 lda,
    i32* jpvt,
    const f64 rcond,
    const f64 svlmax,
    f64* tau,
    i32* rank,
    f64* sval,
    f64* dwork,
    const i32 ldwork,
    i32* info
);

/**
 * @brief Solve augmented system A*x = b, D*x = 0 in least squares sense.
 *
 * Determines vector x which solves the system of linear equations
 * A*x = b, D*x = 0 in the least squares sense, where A is m-by-n,
 * D is n-by-n diagonal, and b is m-vector. Assumes QR factorization
 * with column pivoting of A is available: A*P = Q*R.
 *
 * Uses Givens rotations to annihilate diagonal matrix D, updating
 * upper triangular matrix R and first n elements of Q'*b.
 *
 * @param[in] cond Condition estimation mode:
 *                 'E': use incremental condition estimation
 *                 'N': check diagonal entries for zero
 *                 'U': use rank already stored in RANK
 * @param[in] n Order of matrix R (n >= 0)
 * @param[in,out] r DOUBLE PRECISION array, dimension (ldr,n)
 *                  In: upper triangular matrix R
 *                  Out: full upper triangle unaltered,
 *                       strict lower triangle contains strict upper triangle
 *                       (transposed) of upper triangular matrix S
 * @param[in] ldr Leading dimension of r (ldr >= max(1,n))
 * @param[in] ipvt INTEGER array, dimension (n)
 *                 Permutation matrix P: column j of P is column ipvt[j]
 *                 of identity matrix (1-based indices)
 * @param[in] diag DOUBLE PRECISION array, dimension (n)
 *                 Diagonal elements of matrix D
 * @param[in] qtb DOUBLE PRECISION array, dimension (n)
 *                First n elements of Q'*b
 * @param[in,out] rank INTEGER
 *                     In (COND='U'): numerical rank of S
 *                     Out (COND='E' or 'N'): estimated numerical rank of S
 * @param[out] x DOUBLE PRECISION array, dimension (n)
 *               Least squares solution of A*x = b, D*x = 0
 * @param[in] tol DOUBLE PRECISION
 *                Tolerance for rank determination (COND='E' only)
 *                tol > 0: lower bound for reciprocal condition number
 *                tol <= 0: use default n*eps
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   On exit: first n elements contain diagonal of S,
 *                           next n elements contain solution z
 * @param[in] ldwork Length of dwork
 *                   COND='E': ldwork >= 4*n
 *                   COND!='E': ldwork >= 2*n
 * @param[out] info Exit code (0=success, <0=invalid parameter -info)
 */
void mb02yd(
    const char* cond,
    const i32 n,
    f64* r,
    const i32 ldr,
    const i32* ipvt,
    const f64* diag,
    const f64* qtb,
    i32* rank,
    f64* x,
    const f64 tol,
    f64* dwork,
    const i32 ldwork,
    i32* info
);

/**
 * @brief Minimum norm least squares solution using SVD.
 *
 * Computes the minimum norm least squares solution of one of:
 *   op(R)*X = alpha*B     (SIDE='L')
 *   X*op(R) = alpha*B     (SIDE='R')
 *
 * where alpha is a real scalar, op(R) is R or R', R is an L-by-L real
 * upper triangular matrix, B is an M-by-N real matrix, and L = M for
 * SIDE='L' or L = N for SIDE='R'. Uses singular value decomposition,
 * R = Q*S*P', assuming R may be rank deficient.
 *
 * @param[in] fact 'N' to compute SVD, 'F' if SVD already factored
 * @param[in] side 'L' for op(R)*X = alpha*B, 'R' for X*op(R) = alpha*B
 * @param[in] trans 'N' for op(R)=R, 'T'/'C' for op(R)=R'
 * @param[in] jobp 'P' to compute/use pinv(R), 'N' otherwise
 * @param[in] m Number of rows of B (m >= 0)
 * @param[in] n Number of columns of B (n >= 0)
 * @param[in] alpha Scalar multiplier (if 0, X set to zero)
 * @param[in] rcond Rank threshold: SV(i) <= RCOND*SV(1) treated as zero
 *                  If RCOND <= 0, uses machine epsilon. Not used if FACT='F'.
 * @param[in,out] rank Rank of R. Input if FACT='F', output if FACT='N'.
 * @param[in,out] r DOUBLE PRECISION array, dimension (ldr,L)
 *                  If FACT='N': In: upper triangular R
 *                              Out: P' (scaled by inv(S) if JOBP='P')
 *                  If FACT='F': In: P' from prior SVD
 * @param[in] ldr Leading dimension of r (ldr >= max(1,L))
 * @param[in,out] q DOUBLE PRECISION array, dimension (ldq,L)
 *                  If FACT='N': Out: orthogonal matrix Q
 *                  If FACT='F': In: orthogonal matrix Q from prior SVD
 * @param[in] ldq Leading dimension of q (ldq >= max(1,L))
 * @param[in,out] sv DOUBLE PRECISION array, dimension (L)
 *                   If FACT='N': Out: first RANK entries are 1/SV(1:RANK),
 *                                rest are SV(RANK+1:L) in descending order
 *                   If FACT='F': In: reciprocal SVs and remaining SVs
 * @param[in,out] b DOUBLE PRECISION array, dimension (ldb,n)
 *                  In: matrix B (if alpha != 0)
 *                  Out: solution X
 * @param[in] ldb Leading dimension of b (ldb >= max(1,m))
 * @param[out] rp DOUBLE PRECISION array, dimension (ldrp,L)
 *                If JOBP='P' and RANK > 0: pinv(R)
 *                Not referenced if JOBP='N'
 * @param[in] ldrp Leading dimension of rp (ldrp >= L if JOBP='P', else >= 1)
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   On exit: dwork[0] = optimal ldwork
 * @param[in] ldwork Workspace size (>= max(1,L) if FACT='F', >= max(1,5*L) if FACT='N')
 *                   Optimal: max(1,L,M*N) if FACT='F', max(1,5*L,M*N) if FACT='N'
 *                   If ldwork=-1, workspace query
 * @param[out] info Exit code (0=success, <0=invalid parameter,
 *                  >0=SVD did not converge)
 */
void mb02ud(
    const char* fact, const char* side, const char* trans, const char* jobp,
    const i32 m, const i32 n, const f64 alpha, const f64 rcond,
    i32* rank, f64* r, const i32 ldr, f64* q, const i32 ldq,
    f64* sv, f64* b, const i32 ldb, f64* rp, const i32 ldrp,
    f64* dwork, const i32 ldwork, i32* info
);

/**
 * @brief Matrix rank determination by incremental condition estimation during QR factorization.
 *
 * Computes a rank-revealing QR factorization of a real general M-by-N matrix A,
 * which may be rank-deficient, and estimates its effective rank using incremental
 * condition estimation.
 *
 * The routine uses a truncated QR factorization with column pivoting
 * A * P = Q * R, where R = [R11 R12; 0 R22], with R11 defined as the largest
 * leading upper triangular submatrix whose estimated condition number is less
 * than 1/RCOND. The order of R11, RANK, is the effective rank of A.
 *
 * @param[in] m Number of rows of matrix A. m >= 0.
 * @param[in] n Number of columns of matrix A. n >= 0.
 * @param[in,out] a DOUBLE PRECISION array, dimension (lda, n)
 *                  On entry: M-by-N matrix A
 *                  On exit: RANK-by-RANK upper triangular R11 and QR factorization data
 * @param[in] lda Leading dimension of array A. lda >= max(1,m).
 * @param[in] rcond Threshold for rank determination (0 <= rcond <= 1)
 * @param[in] svlmax Estimate of largest singular value of parent matrix (>= 0)
 * @param[out] rank Effective rank of A
 * @param[out] sval DOUBLE PRECISION array, dimension (3)
 *                  Singular value estimates: [largest, rank-th, (rank+1)-th]
 * @param[out] jpvt INTEGER array, dimension (n)
 *                  Pivot indices (1-based, Fortran style)
 * @param[out] tau DOUBLE PRECISION array, dimension (min(m,n))
 *                 Scalar factors of elementary reflectors
 * @param[out] dwork DOUBLE PRECISION array, dimension (3*n-1)
 * @param[out] info Exit code (0 = success, <0 = invalid parameter)
 */
void mb03oy(i32 m, i32 n, f64* a, i32 lda, f64 rcond,
            f64 svlmax, i32* rank, f64* sval, i32* jpvt,
            f64* tau, f64* dwork, i32* info);

/**
 * @brief Performs a QR factorization update.
 *
 * MB04OW performs the QR factorization
 *
 *      ( U  ) = Q*( R ),  where  U = ( U1  U2 ),  R = ( R1  R2 ),
 *      ( x' )     ( 0 )              ( 0   T  )       ( 0   R3 )
 *
 * where U and R are (m+n)-by-(m+n) upper triangular matrices, x is
 * an m+n element vector, U1 is m-by-m, T is n-by-n, stored
 * separately, and Q is an (m+n+1)-by-(m+n+1) orthogonal matrix.
 *
 * The transformations performed are also applied to the (m+n+1)-by-p
 * matrix ( B' C' d )' (' denotes transposition), where B, C, and d'
 * are m-by-p, n-by-p, and 1-by-p matrices, respectively.
 *
 * @param[in] m The number of rows of the matrix ( U1  U2 ). M >= 0.
 * @param[in] n The order of the matrix T. N >= 0.
 * @param[in] p The number of columns of the matrices B and C. P >= 0.
 * @param[in,out] a Array of dimension (LDA, N+M). On entry, the leading M-by-(M+N) 
 *                  upper trapezoidal part contains ( U1 U2 ). On exit, ( R1 R2 ).
 * @param[in] lda The leading dimension of the array A. LDA >= max(1,M).
 * @param[in,out] t Array of dimension (LDT, N). On entry, the leading N-by-N 
 *                  upper triangular part contains T. On exit, R3.
 * @param[in] ldt The leading dimension of the array T. LDT >= max(1,N).
 * @param[in,out] x Array of dimension (1+(M+N-1)*INCX). On entry, the vector x. 
 *                  On exit, the content is changed (destroyed).
 * @param[in] incx The increment for the elements of X. INCX > 0.
 * @param[in,out] b Array of dimension (LDB, P). On entry, B. On exit, transformed B.
 * @param[in] ldb The leading dimension of the array B. LDB >= max(1,M) if P > 0.
 * @param[in,out] c Array of dimension (LDC, P). On entry, C. On exit, transformed C.
 * @param[in] ldc The leading dimension of the array C. LDC >= max(1,N) if P > 0.
 * @param[in,out] d Array of dimension (1+(P-1)*INCD). On entry, the vector d. 
 *                  On exit, transformed d.
 * @param[in] incd The increment for the elements of D. INCD > 0.
 */
void mb04ow(i32 m, i32 n, i32 p, f64 *a, i32 lda, f64 *t, i32 ldt, 
            f64 *x, i32 incx, f64 *b, i32 ldb, f64 *c, i32 ldc, 
            f64 *d, i32 incd);

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
            f64 *y, i32 ldy, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Compute the Jacobian of the error function for a neural network.
 *
 * Computes the Jacobian of the error function for a neural network of the structure:
 *
 *          - tanh(w1*z+b1) -
 *        /      :            \
 *      z ---    :          --- sum(ws(i)*...)+ b(n+1)  --- y,
 *        \      :            /
 *          - tanh(wn*z+bn) -
 *
 * for the single-output case.
 *
 * @param[in] cjte 'C' to compute J'*e, 'N' to skip.
 * @param[in] nsmp Number of training samples.
 * @param[in] nz Length of each input sample.
 * @param[in] l Length of each output sample (must be 1).
 * @param[in,out] ipar Integer parameters. ipar[0] contains NN.
 * @param[in] lipar Length of ipar.
 * @param[in] wb Weights and biases vector.
 * @param[in] lwb Length of wb.
 * @param[in] z Input samples matrix (nsmp x nz).
 * @param[in] ldz Leading dimension of z.
 * @param[in] e Error vector (nsmp).
 * @param[out] j Jacobian matrix (nsmp x nwb).
 * @param[in] ldj Leading dimension of j.
 * @param[out] jte Vector J'*e (nwb).
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01by(const char *cjte, i32 nsmp, i32 nz, i32 l, i32 *ipar, i32 lipar, 
            const f64 *wb, i32 lwb, const f64 *z, i32 ldz, const f64 *e, 
            f64 *j, i32 ldj, f64 *jte, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief QR factorization with column pivoting for Levenberg-Marquardt.
 *
 * This routine is an interface to SLICOT Library routine MD03BX.
 *
 * @param[in] n Number of columns of Jacobian matrix J.
 * @param[in] ipar Integer parameters. ipar[0] must contain M (rows of J).
 * @param[in] lipar Length of ipar.
 * @param[in] fnorm Euclidean norm of error vector e.
 * @param[in,out] j Jacobian matrix (M x N). On exit, upper triangular R.
 * @param[in,out] ldj Leading dimension of J.
 * @param[in,out] e Error vector (M). On exit, Q'*e.
 * @param[out] jnorms Euclidean norms of columns of J.
 * @param[out] gnorm 1-norm of scaled gradient.
 * @param[out] ipvt Permutation matrix P.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void md03ba(i32 n, const i32 *ipar, i32 lipar, f64 fnorm, f64 *j, i32 *ldj, 
            f64 *e, f64 *jnorms, f64 *gnorm, i32 *ipvt, f64 *dwork, 
            i32 ldwork, i32 *info);

/**
 * @brief Solve system of linear equations R*x = b or R'*x = b in least squares sense.
 *
 * Solves R*x = b or R'*x = b where R is an n-by-n block upper triangular matrix.
 *
 * @param[in] cond Condition estimation mode ('E', 'N', 'U').
 * @param[in] uplo Storage scheme ('U', 'L').
 * @param[in] trans Form of system ('N', 'T', 'C').
 * @param[in] n Order of matrix R.
 * @param[in] ipar Integer parameters (st, bn, bsm, bsn).
 * @param[in] lipar Length of ipar.
 * @param[in,out] r Matrix R (ldr x nc).
 * @param[in] ldr Leading dimension of R.
 * @param[in] sdiag Diagonal elements of blocks (if uplo='L').
 * @param[in] s Transpose of last block column (if uplo='L').
 * @param[in] lds Leading dimension of S.
 * @param[in,out] b Right hand side vector b. On exit, solution x.
 * @param[in,out] ranks Numerical ranks of submatrices.
 * @param[in] tol Tolerance for rank determination.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01br(const char *cond, const char *uplo, const char *trans, i32 n, 
            const i32 *ipar, i32 lipar, f64 *r, i32 ldr, f64 *sdiag, 
            f64 *s, i32 lds, f64 *b, i32 *ranks, f64 tol, f64 *dwork, 
            i32 ldwork, i32 *info);

/**
 * @brief Compute Levenberg-Marquardt parameter for compressed Jacobian.
 *
 * This routine is an interface to SLICOT Library routine MD03BY.
 *
 * @param[in] cond Condition estimation mode ('E', 'N', 'U').
 * @param[in] n Order of matrix R.
 * @param[in] ipar Integer parameters (unused, for compatibility).
 * @param[in] lipar Length of ipar.
 * @param[in,out] r Upper triangular matrix R.
 * @param[in] ldr Leading dimension of R.
 * @param[in] ipvt Permutation matrix P.
 * @param[in] diag Diagonal scaling matrix D.
 * @param[in] qtb First n elements of Q'*b.
 * @param[in] delta Trust region radius.
 * @param[in,out] par Levenberg-Marquardt parameter.
 * @param[in,out] ranks Numerical rank.
 * @param[out] x Least squares solution.
 * @param[out] rx Residual -R*P'*x.
 * @param[in] tol Tolerance for rank estimation.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void md03bb(const char *cond, i32 n, const i32 *ipar, i32 lipar, f64 *r, 
            i32 ldr, const i32 *ipvt, const f64 *diag, const f64 *qtb, 
            f64 delta, f64 *par, i32 *ranks, f64 *x, f64 *rx, f64 tol, 
            f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief QR factorization of Jacobian in compressed form.
 *
 * Computes QR factorization with column pivoting of Jacobian J in compressed form.
 *
 * @param[in] n Number of columns of J.
 * @param[in] ipar Integer parameters (st, bn, bsm, bsn).
 * @param[in] lipar Length of ipar.
 * @param[in] fnorm Norm of error vector.
 * @param[in,out] j Jacobian matrix (ldj x nc).
 * @param[in] ldj Leading dimension of J.
 * @param[in,out] e Error vector.
 * @param[out] jnorms Euclidean norms of columns of J.
 * @param[out] gnorm 1-norm of scaled gradient.
 * @param[out] ipvt Permutation matrix P.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01bs(i32 n, const i32 *ipar, i32 lipar, f64 fnorm, f64 *j, i32 *ldj, 
            f64 *e, f64 *jnorms, f64 *gnorm, i32 *ipvt, f64 *dwork, 
            i32 ldwork, i32 *info);

/**
 * @brief Compute complex Givens rotation in real arithmetic.
 *
 * Computes parameters for complex Givens rotation such that:
 *
 *     (    C      SR+SI*I )   ( XR+XI*I )   ( ZR+ZI*I )
 *     (                   ) * (         ) = (         )
 *     ( -SR+SI*I     C    )   ( YR+YI*I )   (    0    )
 *
 * where C**2 + |SR+SI*I|**2 = 1.
 *
 * Adapted from LAPACK ZLARTG for real data representation.
 * Avoids unnecessary overflow/underflow.
 *
 * @param[in] xr Real part of X
 * @param[in] xi Imaginary part of X
 * @param[in] yr Real part of Y
 * @param[in] yi Imaginary part of Y
 * @param[out] c Cosine parameter (real)
 * @param[out] sr Real part of sine parameter
 * @param[out] si Imaginary part of sine parameter
 * @param[out] zr Real part of result Z
 * @param[out] zi Imaginary part of result Z
 */
void sg03br(
    const f64 xr, const f64 xi, const f64 yr, const f64 yi,
    f64* c, f64* sr, f64* si, f64* zr, f64* zi
);

/**
 * @brief Solve generalized Lyapunov equation for Cholesky factor.
 *
 * Computes Cholesky factor U (X = op(U)^T * op(U)) solving the generalized
 * c-stable continuous-time Lyapunov equation:
 *
 *   op(A)^T * X * op(E) + op(E)^T * X * op(A) = -SCALE^2 * op(B)^T * op(B),
 *
 * or the generalized d-stable discrete-time Lyapunov equation:
 *
 *   op(A)^T * X * op(A) - op(E)^T * X * op(E) = -SCALE^2 * op(B)^T * op(B),
 *
 * where op(K) is either K or K^T. A, E are N-by-N, op(B) is M-by-N.
 * Result U is N-by-N upper triangular with non-negative diagonal entries.
 * SCALE is set to avoid overflow in U.
 *
 * If FACT='N', pencil A-lambda*E is reduced to generalized Schur form.
 * If FACT='F', generalized Schur factors must be supplied on entry.
 *
 * @param[in] dico 'C' for continuous-time, 'D' for discrete-time
 * @param[in] fact 'N' to compute factorization, 'F' if factorization supplied
 * @param[in] trans 'N' for op(K)=K, 'T' for op(K)=K^T
 * @param[in] n Order of matrices A and E (n >= 0)
 * @param[in] m Number of rows in op(B) (m >= 0)
 * @param[in,out] a DOUBLE PRECISION array, dimension (lda,n)
 *                  In: Matrix A (if FACT='F': generalized Schur factor)
 *                  Out: Generalized Schur factor (if FACT='N')
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in,out] e DOUBLE PRECISION array, dimension (lde,n)
 *                  In: Matrix E (if FACT='F': generalized Schur factor)
 *                  Out: Generalized Schur factor (if FACT='N')
 * @param[in] lde Leading dimension of E (lde >= max(1,n))
 * @param[in,out] q DOUBLE PRECISION array, dimension (ldq,n)
 *                  In: Orthogonal matrix Q (if FACT='F')
 *                  Out: Orthogonal matrix Q from factorization (if FACT='N')
 * @param[in] ldq Leading dimension of Q (ldq >= max(1,n))
 * @param[in,out] z DOUBLE PRECISION array, dimension (ldz,n)
 *                  In: Orthogonal matrix Z (if FACT='F')
 *                  Out: Orthogonal matrix Z from factorization (if FACT='N')
 * @param[in] ldz Leading dimension of Z (ldz >= max(1,n))
 * @param[in,out] b DOUBLE PRECISION array, dimension (ldb,n1)
 *                  In: Matrix B (size depends on TRANS)
 *                  Out: Cholesky factor U
 * @param[in] ldb Leading dimension of B
 * @param[out] scale Scaling factor (0 < scale <= 1)
 * @param[out] alphar DOUBLE PRECISION array, dimension (n)
 *                    Real parts of eigenvalues of pencil A-lambda*E
 * @param[out] alphai DOUBLE PRECISION array, dimension (n)
 *                    Imaginary parts of eigenvalues of pencil A-lambda*E
 * @param[out] beta DOUBLE PRECISION array, dimension (n)
 *                  Scaling factors for eigenvalues
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   Workspace, dwork[0] returns optimal ldwork
 * @param[in] ldwork Workspace size (ldwork >= max(1,4*n,6*n-6) if FACT='N',
 *                   ldwork >= max(1,2*n,6*n-6) if FACT='F')
 *                   If ldwork=-1, workspace query
 * @param[out] info Exit code (0=success, <0=invalid parameter, 1=singular,
 *                  2=not quasitriangular, 3=eigenvalues not conjugate,
 *                  4=factorization failed, 5=not c-stable, 6=not d-stable,
 *                  7=DSYEVX failed)
 */
void sg03bd(
    const char* dico, const char* fact, const char* trans,
    const i32 n, const i32 m,
    f64* a, const i32 lda,
    f64* e, const i32 lde,
    f64* q, const i32 ldq,
    f64* z, const i32 ldz,
    f64* b, const i32 ldb,
    f64* scale,
    f64* alphar, f64* alphai, f64* beta,
    f64* dwork, const i32 ldwork,
    i32* info
);

/**
 * @brief Solve generalized Sylvester equation for small systems.
 *
 * Solves for X the generalized Sylvester equation:
 *
 *     A^T * X * C + E^T * X * D = SCALE * Y,    (TRANS='N')
 *
 * or the transposed equation:
 *
 *     A * X * C^T + E * X * D^T = SCALE * Y,    (TRANS='T')
 *
 * where A and E are M-by-M matrices (A upper quasitriangular, E upper triangular),
 * C and D are N-by-N matrices, X and Y are M-by-N matrices. N must be 1 or 2.
 * The pencil A - lambda*E must be in generalized real Schur form.
 * SCALE is set to avoid overflow in X.
 *
 * @param[in] trans 'N' for equation (1), 'T' for transposed equation
 * @param[in] m Order of matrices A and E (m >= 0)
 * @param[in] n Order of matrices C and D (n = 1 or 2)
 * @param[in] a DOUBLE PRECISION array, dimension (lda,m)
 *              Upper quasitriangular matrix A
 * @param[in] lda Leading dimension of A (lda >= max(1,m))
 * @param[in] c DOUBLE PRECISION array, dimension (ldc,n)
 *              Matrix C
 * @param[in] ldc Leading dimension of C (ldc >= max(1,n))
 * @param[in] e DOUBLE PRECISION array, dimension (lde,m)
 *              Upper triangular matrix E
 * @param[in] lde Leading dimension of E (lde >= max(1,m))
 * @param[in] d DOUBLE PRECISION array, dimension (ldd,n)
 *              Matrix D
 * @param[in] ldd Leading dimension of D (ldd >= max(1,n))
 * @param[in,out] x DOUBLE PRECISION array, dimension (ldx,n)
 *                  In: Right-hand side Y
 *                  Out: Solution X
 * @param[in] ldx Leading dimension of X (ldx >= max(1,m))
 * @param[out] scale Scaling factor (0 < scale <= 1)
 * @param[out] info Exit code (0 = success, <0 = invalid parameter,
 *                  1 = nearly singular, perturbed values used)
 */
void sg03bw(
    const char* trans,
    const i32 m, const i32 n,
    const f64* a, const i32 lda,
    const f64* c, const i32 ldc,
    const f64* e, const i32 lde,
    const f64* d, const i32 ldd,
    f64* x, const i32 ldx,
    f64* scale,
    i32* info
);

/**
 * @brief Solve generalized discrete-time Lyapunov equation for Cholesky factor.
 *
 * Computes Cholesky factor U (X = U^T * U or X = U * U^T) solving
 * generalized d-stable discrete-time Lyapunov equation:
 *
 *   TRANS='N': A^T * X * A - E^T * X * E = -SCALE^2 * B^T * B
 *   TRANS='T': A * X * A^T - E * X * E^T = -SCALE^2 * B * B^T
 *
 * A quasitriangular, E and B upper triangular. Pencil A-lambda*E must be
 * d-stable (eigenvalue moduli < 1).
 *
 * @param[in] trans 'N' for equation (1), 'T' for equation (2)
 * @param[in] n Order of matrices (N >= 0)
 * @param[in] a DOUBLE PRECISION array, dimension (lda,N), quasitriangular A
 * @param[in] lda Leading dimension of A (lda >= max(1,N))
 * @param[in] e DOUBLE PRECISION array, dimension (lde,N), upper triangular E
 * @param[in] lde Leading dimension of E (lde >= max(1,N))
 * @param[in,out] b DOUBLE PRECISION array, dimension (ldb,N)
 *                  On entry: upper triangular B
 *                  On exit: Cholesky factor U
 * @param[in] ldb Leading dimension of B (ldb >= max(1,N))
 * @param[out] scale Scaling factor (0 < scale <= 1)
 * @param[out] dwork DOUBLE PRECISION workspace, dimension (6*N-6)
 * @param[out] info Exit code (0=success, <0=parameter error, 1=near singular,
 *                  2=not complex conjugate, 3=not d-stable, 4=DSYEVX failed)
 */
void sg03bu(
    const char* trans, const i32 n,
    const f64* a, const i32 lda,
    const f64* e, const i32 lde,
    f64* b, const i32 ldb,
    f64* scale, f64* dwork,
    i32* info
);

/**
 * @brief Solve generalized continuous-time Lyapunov equation for Cholesky factor.
 *
 * Computes Cholesky factor U (X = U^T * U or X = U * U^T) solving
 * generalized c-stable continuous-time Lyapunov equation:
 *
 *   TRANS='N': A^T * X * E + E^T * X * A = -SCALE^2 * B^T * B
 *   TRANS='T': A * X * E^T + E * X * A^T = -SCALE^2 * B * B^T
 *
 * A quasitriangular, E and B upper triangular. Pencil A-lambda*E must be
 * c-stable (eigenvalues with negative real parts).
 *
 * @param[in] trans 'N' for equation (1), 'T' for equation (2)
 * @param[in] n Order of matrices (N >= 0)
 * @param[in] a DOUBLE PRECISION array, dimension (lda,N), quasitriangular A
 * @param[in] lda Leading dimension of A (lda >= max(1,N))
 * @param[in] e DOUBLE PRECISION array, dimension (lde,N), upper triangular E
 * @param[in] lde Leading dimension of E (lde >= max(1,N))
 * @param[in,out] b DOUBLE PRECISION array, dimension (ldb,N)
 *                  On entry: upper triangular B
 *                  On exit: Cholesky factor U
 * @param[in] ldb Leading dimension of B (ldb >= max(1,N))
 * @param[out] scale Scaling factor (0 < scale <= 1)
 * @param[out] dwork DOUBLE PRECISION workspace, dimension (6*N-6)
 * @param[out] info Exit code (0=success, <0=parameter error, 1=near singular,
 *                  2=not complex conjugate, 3=not c-stable)
 */
void sg03bv(
    const char* trans, const i32 n,
    const f64* a, const i32 lda,
    const f64* e, const i32 lde,
    f64* b, const i32 ldb,
    f64* scale, f64* dwork,
    i32* info
);

/**
 * @brief Solve 2-by-2 generalized Lyapunov equation.
 *
 * Solves for Cholesky factor U (X = op(U)^T * op(U)) the generalized
 * continuous-time or discrete-time Lyapunov equation:
 *
 *   Continuous (DICO='C'):
 *     op(A)^T * X * op(E) + op(E)^T * X * op(A) = -SCALE^2 * op(B)^T * op(B)
 *
 *   Discrete (DICO='D'):
 *     op(A)^T * X * op(A) - op(E)^T * X * op(E) = -SCALE^2 * op(B)^T * op(B)
 *
 * where op(K) = K or K^T, A,B,E,U are 2x2 real matrices, E and B upper triangular.
 * Pencil A-lambda*E must have complex conjugate eigenvalues in the stability region.
 * Also computes auxiliary matrices M1 and M2.
 *
 * @param[in] dico 'C' for continuous-time, 'D' for discrete-time
 * @param[in] trans 'N' for op(K)=K, 'T' for op(K)=K^T
 * @param[in] a DOUBLE PRECISION array, dimension (lda,2), matrix A
 * @param[in] lda Leading dimension of A (lda >= 2)
 * @param[in] e DOUBLE PRECISION array, dimension (lde,2), upper triangular E
 * @param[in] lde Leading dimension of E (lde >= 2)
 * @param[in] b DOUBLE PRECISION array, dimension (ldb,2), upper triangular B
 * @param[in] ldb Leading dimension of B (ldb >= 2)
 * @param[out] u DOUBLE PRECISION array, dimension (ldu,2), Cholesky factor
 * @param[in] ldu Leading dimension of U (ldu >= 2)
 * @param[out] scale Scaling factor (0 < scale <= 1)
 * @param[out] m1 DOUBLE PRECISION array, dimension (ldm1,2), auxiliary matrix
 * @param[in] ldm1 Leading dimension of M1 (ldm1 >= 2)
 * @param[out] m2 DOUBLE PRECISION array, dimension (ldm2,2), auxiliary matrix
 * @param[in] ldm2 Leading dimension of M2 (ldm2 >= 2)
 * @param[out] info Exit code (0=success, 2=not complex conjugate,
 *                  3=eigenvalues not stable, 4=ZSTEIN failed)
 */
void sg03bx(
    const char* dico, const char* trans,
    const f64* a, const i32 lda,
    const f64* e, const i32 lde,
    const f64* b, const i32 ldb,
    f64* u, const i32 ldu,
    f64* scale,
    f64* m1, const i32 ldm1,
    f64* m2, const i32 ldm2,
    i32* info
);

/**
 * @brief Orthogonal reduction of descriptor system to SVD-like coordinate form.
 *
 * Computes orthogonal transformation matrices Q and Z such that the transformed
 * descriptor system (Q'*A*Z - lambda*Q'*E*Z, Q'*B, C*Z) is in SVD-like form:
 *
 *            ( A11  A12 )             ( Er  0 )
 *   Q'*A*Z = (          ) ,  Q'*E*Z = (       )
 *            ( A21  A22 )             (  0  0 )
 *
 * where Er is upper triangular and invertible. Optionally reduces A22 to:
 *
 *        ( Ar  X )                ( Ar  0 )
 *  A22 = (       )  (JOBA='T') or (       )  (JOBA='R')
 *        (  0  0 )                (  0  0 )
 *
 * with Ar upper triangular invertible, X full or zero.
 *
 * @param[in] compq Controls Q computation:
 *                  'N' = do not compute Q
 *                  'I' = initialize Q to identity and return Q
 *                  'U' = update existing Q (Q := Q1*Q)
 * @param[in] compz Controls Z computation:
 *                  'N' = do not compute Z
 *                  'I' = initialize Z to identity and return Z
 *                  'U' = update existing Z (Z := Z1*Z)
 * @param[in] joba Controls A22 reduction:
 *                 'N' = do not reduce A22
 *                 'R' = reduce A22 to SVD-like upper triangular form
 *                 'T' = reduce A22 to upper trapezoidal form
 * @param[in] l Number of rows of A, B, E (L >= 0)
 * @param[in] n Number of columns of A, E, C (N >= 0)
 * @param[in] m Number of columns of B (M >= 0)
 * @param[in] p Number of rows of C (P >= 0)
 * @param[in,out] a DOUBLE PRECISION array, dimension (lda,n)
 *                  In: L-by-N state dynamics matrix A
 *                  Out: Transformed matrix Q'*A*Z
 * @param[in] lda Leading dimension of A (lda >= max(1,l))
 * @param[in,out] e DOUBLE PRECISION array, dimension (lde,n)
 *                  In: L-by-N descriptor matrix E
 *                  Out: Transformed matrix Q'*E*Z in SVD-like form
 * @param[in] lde Leading dimension of E (lde >= max(1,l))
 * @param[in,out] b DOUBLE PRECISION array, dimension (ldb,m)
 *                  In: L-by-M input/state matrix B
 *                  Out: Transformed matrix Q'*B
 * @param[in] ldb Leading dimension of B (ldb >= max(1,l) if m>0, else >= 1)
 * @param[in,out] c DOUBLE PRECISION array, dimension (ldc,n)
 *                  In: P-by-N state/output matrix C
 *                  Out: Transformed matrix C*Z
 * @param[in] ldc Leading dimension of C (ldc >= max(1,p))
 * @param[in,out] q DOUBLE PRECISION array, dimension (ldq,l)
 *                  If compq='I': Out: orthogonal matrix Q
 *                  If compq='U': In: Q1, Out: Q1*Q
 *                  If compq='N': Not referenced
 * @param[in] ldq Leading dimension of Q (ldq >= max(1,l) if compq!='N', else >= 1)
 * @param[in,out] z DOUBLE PRECISION array, dimension (ldz,n)
 *                  If compz='I': Out: orthogonal matrix Z
 *                  If compz='U': In: Z1, Out: Z1*Z
 *                  If compz='N': Not referenced
 * @param[in] ldz Leading dimension of Z (ldz >= max(1,n) if compz!='N', else >= 1)
 * @param[out] ranke Rank of matrix E (order of Er)
 * @param[out] rnka22 Rank of A22 (order of Ar, if joba='R' or 'T')
 * @param[in] tol Tolerance for rank determination (0 < tol < 1)
 *                If tol <= 0, uses default: l*n*eps
 * @param[out] iwork INTEGER array, dimension (n)
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   On exit, dwork[0] returns optimal ldwork
 * @param[in] ldwork Length of dwork (>= max(1, n+p, min(l,n)+max(3*n-1,m,l)))
 *                   If ldwork=-1, workspace query (returns optimal size in dwork[0])
 * @param[out] info Exit code (0 = success, <0 = invalid parameter -i)
 */
void tg01fd(
    const char* compq, const char* compz, const char* joba,
    const i32 l, const i32 n, const i32 m, const i32 p,
    f64* a, const i32 lda,
    f64* e, const i32 lde,
    f64* b, const i32 ldb,
    f64* c, const i32 ldc,
    f64* q, const i32 ldq,
    f64* z, const i32 ldz,
    i32* ranke, i32* rnka22,
    const f64 tol,
    i32* iwork, f64* dwork, const i32 ldwork,
    i32* info
);

/**
 * @brief Transpose all or part of a matrix.
 *
 * Transposes an M-by-N matrix A into N-by-M matrix B.
 * Supports full, upper triangular, or lower triangular transposition.
 *
 * @param[in] job Specifies part to transpose:
 *                'U' = upper triangular/trapezoidal part only
 *                'L' = lower triangular/trapezoidal part only
 *                Otherwise = full matrix
 * @param[in] m Number of rows of A (m >= 0)
 * @param[in] n Number of columns of A (n >= 0)
 * @param[in] a Input matrix, dimension (lda,n), column-major
 * @param[in] lda Leading dimension of A (lda >= max(1,m))
 * @param[out] b Output matrix (transpose), dimension (ldb,m), column-major
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n))
 */
void ma02ad(const char* job, const i32 m, const i32 n,
            const f64* a, const i32 lda,
            f64* b, const i32 ldb);

/**
 * @brief Store by symmetry the upper or lower triangle of a symmetric matrix.
 *
 * Completes a symmetric matrix by copying one triangle to the other.
 * Given upper triangle, constructs lower triangle (or vice versa).
 *
 * @param[in] uplo Specifies which part is given:
 *                 'U' = upper triangular part given
 *                 'L' = lower triangular part given
 *                 Other values = no operation
 * @param[in] n Order of matrix A (n >= 0)
 * @param[in,out] a Matrix, dimension (lda,n), column-major
 *                  In: N-by-N upper/lower triangle contains data
 *                  Out: Full N-by-N symmetric matrix
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 */
void ma02ed(const char uplo, i32 n, f64 *a, i32 lda);

/**
 * @brief Compute coefficients for modified hyperbolic plane rotation.
 *
 * Computes c and s (c^2 + s^2 = 1) such that:
 *     y1 = (1/c) * x1 - (s/c) * x2 = sqrt(x1^2 - x2^2)
 *     y2 = -s * y1 + c * x2 = 0
 *
 * The input must satisfy either x1 = x2 = 0, or abs(x2) < abs(x1).
 *
 * @param[in,out] x1 On entry: x1. On exit: y1 = sqrt(x1^2 - x2^2)
 * @param[in] x2 The value x2
 * @param[out] c Cosine of the modified hyperbolic rotation
 * @param[out] s Sine of the modified hyperbolic rotation
 * @param[out] info 0 = success, 1 = abs(x2) >= abs(x1) with x1 != 0 or x2 != 0
 */
void ma02fd(f64 *x1, f64 x2, f64 *c, f64 *s, i32 *info);

/**
 * @brief QR factorization of special structured block matrix.
 *
 * Computes QR factorization of first block column and applies orthogonal
 * transformations to second block column:
 *         [[R], [A B]] -> Q' * [[R], [A B]] = [[R_bar C], [0 D]]
 * where R and R_bar are upper triangular.
 *
 * @param[in] uplo Indicates structure of A:
 *                 'U' = upper trapezoidal/triangular
 *                 'F' = full matrix
 * @param[in] n Order of matrices R and R_bar (n >= 0)
 * @param[in] m Number of columns of B, C, D (m >= 0)
 * @param[in] p Number of rows of A, B, D (p >= 0)
 * @param[in,out] r Matrix R, dimension (ldr,n)
 *                  In: n-by-n upper triangular R
 *                  Out: n-by-n upper triangular R_bar
 * @param[in] ldr Leading dimension of r (ldr >= max(1,n))
 * @param[in,out] a Matrix A, dimension (lda,n)
 *                  In: p-by-n matrix (full or upper trapezoidal)
 *                  Out: Householder vectors v_i
 * @param[in] lda Leading dimension of a (lda >= max(1,p))
 * @param[in,out] b Matrix B, dimension (ldb,m)
 *                  In: p-by-m matrix B
 *                  Out: p-by-m matrix D
 * @param[in] ldb Leading dimension of b (ldb >= max(1,p))
 * @param[out] c Matrix C, dimension (ldc,m)
 *               n-by-m matrix C
 * @param[in] ldc Leading dimension of c (ldc >= max(1,n))
 * @param[out] tau Scalar factors of elementary reflectors, dimension (n)
 * @param[out] dwork Workspace, dimension (n)
 */
void mb04kd(const char uplo, i32 n, i32 m, i32 p,
            f64 *r, i32 ldr,
            f64 *a, i32 lda,
            f64 *b, i32 ldb,
            f64 *c, i32 ldc,
            f64 *tau,
            f64 *dwork);

/**
 * @brief QR factorization of matrix with lower-left zero triangle
 *
 * Computes QR factorization A = Q*R of n-by-m matrix A having p-by-min(p,m)
 * zero triangle in lower left corner. Optionally applies Q' to n-by-l matrix B.
 * Exploits structure for efficiency (useful in Kalman filtering).
 *
 * Example structure (n=8, m=7, p=2):
 *     [ x x x x x x x ]
 *     [ x x x x x x x ]
 *     [ x x x x x x x ]
 *     [ x x x x x x x ]
 * A = [ x x x x x x x ]
 *     [ x x x x x x x ]
 *     [ 0 x x x x x x ]  <- p rows with
 *     [ 0 0 x x x x x ]  <- lower-left zeros
 *
 * @param[in] n Number of rows of A (n >= 0)
 * @param[in] m Number of columns of A (m >= 0)
 * @param[in] p Order of zero triangle (p >= 0)
 * @param[in] l Number of columns of B (l >= 0)
 * @param[in,out] a Matrix A, dimension (lda,m)
 *                  In: n-by-m matrix with p-by-min(p,m) zero lower-left triangle
 *                  Out: min(n,m)-by-m upper trapezoidal R, Householder vectors below diagonal
 * @param[in] lda Leading dimension of a (lda >= max(1,n))
 * @param[in,out] b Matrix B, dimension (ldb,l)
 *                  In: n-by-l matrix B
 *                  Out: Q'*B
 *                  Not referenced if l=0
 * @param[in] ldb Leading dimension of b (ldb >= max(1,n) if l>0, ldb >= 1 if l=0)
 * @param[out] tau Householder scalar factors, dimension (min(n,m))
 * @param[out] dwork Workspace, dimension (ldwork)
 *                   dwork[0] returns optimal ldwork
 * @param[in] ldwork Workspace size (ldwork >= max(1,m-1,m-p,l))
 *                   If ldwork=-1, workspace query
 * @param[out] info Exit code (0=success, <0=invalid parameter)
 */
void mb04id(i32 n, i32 m, i32 p, i32 l, f64 *a, i32 lda, f64 *b, i32 ldb,
            f64 *tau, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief QR factorization with column pivoting for Levenberg-Marquardt
 *
 * Computes QR factorization with column pivoting of m-by-n matrix J (m >= n):
 * J*P = Q*R, where Q has orthogonal columns, P is permutation, R is upper
 * triangular with diagonal elements of nonincreasing magnitude.
 * Applies Q' to error vector e in-place.
 *
 * @param[in] m Number of rows of Jacobian matrix J (m >= 0)
 * @param[in] n Number of columns of J (m >= n >= 0)
 * @param[in] fnorm Euclidean norm of error vector e (fnorm >= 0)
 * @param[in,out] j Jacobian matrix, dimension (ldj, n)
 *                  In: m-by-n Jacobian matrix
 *                  Out: n-by-n upper triangular R with ldj=n
 * @param[in,out] ldj Leading dimension of J
 *                    In: ldj >= max(1,m)
 *                    Out: ldj = max(1,n)
 * @param[in,out] e Error vector, dimension (m)
 *                  In: Error vector e
 *                  Out: Transformed vector Q'*e
 * @param[out] jnorms Column norms of J (original order), dimension (n)
 * @param[out] gnorm 1-norm of scaled gradient J'*Q'*e/fnorm
 *                   (each element divided by jnorms)
 * @param[out] ipvt Permutation indices, dimension (n)
 *                  Column j of P is column ipvt[j] of identity
 * @param[out] dwork Workspace, dimension (ldwork)
 *                   dwork[0] returns optimal ldwork
 * @param[in] ldwork Workspace size
 *                   ldwork >= 1 if n=0 or m=1
 *                   ldwork >= 4*n+1 if n>1
 * @param[out] info Exit code (0=success, <0=invalid parameter)
 */
void md03bx(
    i32 m, i32 n, f64 fnorm,
    f64* j, i32* ldj, f64* e,
    f64* jnorms, f64* gnorm, i32* ipvt,
    f64* dwork, i32 ldwork, i32* info
);

/**
 * @brief Compute Levenberg-Marquardt parameter for trust region subproblem.
 *
 * Determines parameter PAR such that if x solves the system
 *   A*x = b, sqrt(PAR)*D*x = 0
 * in the least squares sense, then ||D*x|| satisfies the trust region constraint:
 *   either PAR=0 and (||D*x|| - DELTA) <= 0.1*DELTA,
 *   or PAR>0 and abs(||D*x|| - DELTA) <= 0.1*DELTA.
 *
 * Assumes QR factorization A*P = Q*R is available (R, IPVT, Q'*b).
 * Provides upper triangular S such that P'*(A'*A + PAR*D*D)*P = S'*S.
 *
 * @param[in] cond Condition estimation mode:
 *                 'E' = estimate condition of R and S
 *                 'N' = check diagonal entries for zeros only
 *                 'U' = use rank already in RANK parameter
 * @param[in] n Order of matrix R (n >= 0)
 * @param[in,out] r Upper triangular matrix, dimension (ldr,n)
 *                  In: QR factor from A*P=Q*R
 *                  Out: strict lower triangle contains S' (transposed)
 * @param[in] ldr Leading dimension of R (ldr >= max(1,n))
 * @param[in] ipvt Permutation from QR, dimension (n)
 *                 Column j of P is column IPVT(j) of identity
 * @param[in] diag Diagonal scaling matrix D, dimension (n)
 *                 All elements must be nonzero
 * @param[in] qtb First n elements of Q'*b, dimension (n)
 * @param[in] delta Trust region radius (delta > 0)
 * @param[in,out] par Levenberg-Marquardt parameter
 *                    In: initial estimate (par >= 0)
 *                    Out: final estimate
 * @param[in,out] rank Numerical rank
 *                     In: rank of R if COND='U'
 *                     Out: rank of S
 * @param[out] x Least squares solution, dimension (n)
 * @param[out] rx Residual -R*P'*x, dimension (n)
 * @param[in] tol Tolerance for rank estimation if COND='E'
 *                If tol <= 0, use n*eps (machine precision)
 *                Not used if COND='N' or 'U'
 * @param[out] dwork Workspace, dimension (ldwork)
 *                   First n elements contain diagonal of S on exit
 * @param[in] ldwork Workspace size
 *                   ldwork >= 4*n if COND='E'
 *                   ldwork >= 2*n if COND='N' or 'U'
 * @param[out] info Exit code: 0=success, <0=invalid parameter
 */
void md03by(
    const char* cond,
    const i32 n,
    f64* r,
    const i32 ldr,
    const i32* ipvt,
    const f64* diag,
    const f64* qtb,
    const f64 delta,
    f64* par,
    i32* rank,
    f64* x,
    f64* rx,
    const f64 tol,
    f64* dwork,
    const i32 ldwork,
    i32* info
);

/**
 * @brief Levenberg-Marquardt nonlinear least squares optimizer.
 *
 * Minimize sum of squares of m nonlinear functions in n variables
 * using modified Levenberg-Marquardt algorithm with trust region.
 * Requires user-provided FCN (function/Jacobian), QRFACT (QR factorization),
 * and LMPARM (L-M parameter computation) subroutines.
 *
 * @param[in] xinit 'R'=random initialization, 'G'=use given X
 * @param[in] scale 'I'=internal scaling, 'S'=use specified DIAG
 * @param[in] cond 'E'=use condition estimation, 'N'=check diagonal only
 * @param[in] fcn Function pointer for error functions and Jacobian
 * @param[in] qrfact Function pointer for QR factorization with pivoting
 * @param[in] lmparm Function pointer for Levenberg-Marquardt parameter
 * @param[in] m Number of functions (m >= 0)
 * @param[in] n Number of variables (m >= n >= 0)
 * @param[in] itmax Maximum iterations (itmax >= 0)
 * @param[in] factor Initial step bound factor (factor > 0, typically 100)
 * @param[in] nprint Print frequency (IFLAG=0 calls). If <= 0, no printing
 * @param[in] ipar INTEGER array, dimension (lipar). Problem parameters
 * @param[in] lipar Length of ipar (lipar >= 5)
 * @param[in] dpar1 DOUBLE PRECISION array. First parameter set
 * @param[in] ldpar1 Leading dimension/length of dpar1
 * @param[in] dpar2 DOUBLE PRECISION array. Second parameter set
 * @param[in] ldpar2 Leading dimension/length of dpar2
 * @param[in,out] x DOUBLE PRECISION array, dimension (n)
 *                  Input: initial guess (if xinit='G')
 *                  Output: solution minimizing sum of squares
 * @param[in,out] diag DOUBLE PRECISION array, dimension (n)
 *                     Input: scaling factors (if scale='S')
 *                     Output: final scaling factors used
 * @param[out] nfev Number of function evaluations (IFLAG=1)
 * @param[out] njev Number of Jacobian evaluations (IFLAG=2)
 * @param[in] ftol Relative error tolerance for sum of squares
 *                 (ftol < 0 uses sqrt(eps))
 * @param[in] xtol Relative error tolerance for solution
 *                 (xtol < 0 uses sqrt(eps))
 * @param[in] gtol Orthogonality tolerance between e and J columns
 *                 (gtol < 0 uses eps)
 * @param[in] tol Tolerance for rank determination if cond='E'
 *                (tol <= 0 uses n*eps)
 * @param[out] iwork INTEGER array, dimension (n+r)
 *                   iwork[0:n-1]: permutation defining J*P = Q*R
 *                   iwork[n:n+r-1]: ranks of submatrices
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   dwork[0]: optimal ldwork
 *                   dwork[1]: final residual norm
 *                   dwork[2]: iterations performed
 *                   dwork[3]: final Levenberg-Marquardt parameter
 * @param[in] ldwork Length of dwork (see MD03BD documentation)
 * @param[out] iwarn Warning indicator
 *                   <0: user set IFLAG=iwarn
 *                   1: both actual/predicted reductions <= ftol
 *                   2: relative error between iterates <= xtol
 *                   3: conditions 1 and 2 both hold
 *                   4: cosine(e,J) <= gtol
 *                   5: iterations reached itmax
 *                   6: ftol too small
 *                   7: xtol too small
 *                   8: gtol too small
 * @param[out] info Exit code
 *                  0: success
 *                  <0: invalid parameter -info
 *                  1: FCN returned info != 0 for IFLAG=1
 *                  2: FCN returned info != 0 for IFLAG=2
 *                  3: QRFACT returned info != 0
 *                  4: LMPARM returned info != 0
 */
void md03bd(
    const char* xinit,
    const char* scale,
    const char* cond,
    void (*fcn)(i32*, i32, i32, i32*, i32, const f64*, i32, const f64*, i32,
                const f64*, i32*, f64*, f64*, i32*, f64*, i32, i32*),
    void (*qrfact)(i32, const i32*, i32, f64, f64*, i32*, f64*, f64*, f64*,
                   i32*, f64*, i32, i32*),
    void (*lmparm)(const char*, i32, const i32*, i32, f64*, i32, const i32*,
                   const f64*, const f64*, f64, f64*, i32*, f64*, f64*, f64,
                   f64*, i32, i32*),
    i32 m,
    i32 n,
    i32 itmax,
    f64 factor,
    i32 nprint,
    i32* ipar,
    i32 lipar,
    const f64* dpar1,
    i32 ldpar1,
    const f64* dpar2,
    i32 ldpar2,
    f64* x,
    f64* diag,
    i32* nfev,
    i32* njev,
    f64 ftol,
    f64 xtol,
    f64 gtol,
    f64 tol,
    i32* iwork,
    f64* dwork,
    i32 ldwork,
    i32* iwarn,
    i32* info
);

/**
 * @brief Reduce state matrix A to real Schur form via orthogonal transformation
 *
 * Reduces system state matrix A to upper real Schur form by orthogonal similarity
 * transformation A <- U'*A*U, and applies transformation to B <- U'*B and C <- C*U.
 *
 * @param[in] n Order of state matrix A (n >= 0)
 * @param[in] m Number of inputs (columns of B, m >= 0)
 * @param[in] p Number of outputs (rows of C, p >= 0)
 * @param[in,out] a State matrix, dimension (lda,n)
 *                  On entry: original state dynamics matrix
 *                  On exit: real Schur form U'*A*U
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in,out] b Input matrix, dimension (ldb,m)
 *                  On entry: original input matrix
 *                  On exit: transformed matrix U'*B
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n))
 * @param[in,out] c Output matrix, dimension (ldc,n)
 *                  On entry: original output matrix
 *                  On exit: transformed matrix C*U
 * @param[in] ldc Leading dimension of C (ldc >= max(1,p))
 * @param[out] u Orthogonal transformation matrix, dimension (ldu,n)
 *               Schur vectors of A
 * @param[in] ldu Leading dimension of U (ldu >= max(1,n))
 * @param[out] wr Real parts of eigenvalues, dimension (n)
 * @param[out] wi Imaginary parts of eigenvalues, dimension (n)
 * @param[out] dwork Workspace array, dimension (ldwork)
 * @param[in] ldwork Workspace size (ldwork >= 3*n, larger for optimal performance)
 * @param[out] info Exit code: 0=success, <0=invalid parameter, >0=QR failed
 */
void tb01wd(
    const i32 n, const i32 m, const i32 p,
    f64* a, const i32 lda,
    f64* b, const i32 ldb,
    f64* c, const i32 ldc,
    f64* u, const i32 ldu,
    f64* wr, f64* wi,
    f64* dwork, const i32 ldwork,
    i32* info
);

/**
 * @brief Convert discrete-time system to output normal form.
 *
 * Converts a stable discrete-time system (A, B, C, D) with initial state x0
 * into the output normal form, producing parameter vector THETA.
 *
 * The parameter vector THETA contains:
 * - THETA[0:N*L-1]: parameters for A and C matrices
 * - THETA[N*L:N*(L+M)-1]: parameters for B matrix
 * - THETA[N*(L+M):N*(L+M)+L*M-1]: parameters for D matrix
 * - THETA[N*(L+M)+L*M:N*(L+M+1)+L*M-1]: initial state x0
 *
 * Algorithm:
 * 1. Solve Lyapunov equation A'*Q*A - Q = -scale^2*C'*C in Cholesky factor T
 * 2. Transform system using T
 * 3. QR factorization of transposed observability matrix
 * 4. Extract parameters via N orthogonal transformations
 *
 * @param[in] apply Bijective mapping mode:
 *                  'A' = apply bijective mapping to remove norm(THETA_i) < 1 constraint
 *                  'N' = no bijective mapping
 * @param[in] n System order (N >= 0)
 * @param[in] m Number of inputs (M >= 0)
 * @param[in] l Number of outputs (L >= 0)
 * @param[in,out] a State matrix, dimension (LDA,N), column-major.
 *                  On entry: original system matrix (must be stable).
 *                  On exit: transformed system matrix.
 * @param[in] lda Leading dimension of A (>= max(1,N))
 * @param[in,out] b Input matrix, dimension (LDB,M), column-major.
 *                  On entry: original input matrix.
 *                  On exit: transformed input matrix.
 * @param[in] ldb Leading dimension of B (>= max(1,N))
 * @param[in,out] c Output matrix, dimension (LDC,N), column-major.
 *                  On entry: original output matrix.
 *                  On exit: transformed output matrix.
 * @param[in] ldc Leading dimension of C (>= max(1,L))
 * @param[in] d Feedthrough matrix, dimension (LDD,M), column-major (read-only)
 * @param[in] ldd Leading dimension of D (>= max(1,L))
 * @param[in,out] x0 Initial state vector, dimension (N).
 *                   On entry: original initial state.
 *                   On exit: transformed initial state.
 * @param[out] theta Parameter vector, dimension (LTHETA)
 * @param[in] ltheta Length of THETA array (>= N*(L+M+1)+L*M)
 * @param[out] scale Scale factor from Lyapunov equation solver
 * @param[out] dwork Workspace array, dimension (LDWORK)
 * @param[in] ldwork Length of DWORK
 * @param[out] info Exit code:
 *                  = 0: success
 *                  < 0: if INFO = -i, the i-th argument had an illegal value
 *                  = 1: Lyapunov equation could only be solved with scale = 0
 *                  = 2: matrix A is not discrete-time stable
 *                  = 3: QR algorithm failed to converge for matrix A
 */
void tb01vd(const char* apply, i32 n, i32 m, i32 l, f64* a, i32 lda,
            f64* b, i32 ldb, f64* c, i32 ldc, const f64* d, i32 ldd,
            f64* x0, f64* theta, i32 ltheta, f64* scale,
            f64* dwork, i32 ldwork, i32* info);

/**
 * @brief Convert output normal form to state-space representation.
 *
 * Converts a discrete-time system from output normal form (parameter vector THETA)
 * to standard state-space representation (A, B, C, D) with initial state x0.
 *
 * The parameter vector THETA contains:
 * - THETA[0:N*L-1]: parameters for A and C matrices
 * - THETA[N*L:N*(L+M)-1]: parameters for B matrix
 * - THETA[N*(L+M):N*(L+M)+L*M-1]: parameters for D matrix
 * - THETA[N*(L+M)+L*M:N*(L+M+1)+L*M-1]: initial state x0
 *
 * @param[in] apply Bijective mapping mode:
 *                  'A' = apply bijective mapping to remove norm(THETA_i) < 1 constraint
 *                  'N' = no bijective mapping
 * @param[in] n System order (N >= 0)
 * @param[in] m Number of inputs (M >= 0)
 * @param[in] l Number of outputs (L >= 0)
 * @param[in] theta Parameter vector, dimension (LTHETA)
 * @param[in] ltheta Length of THETA array (>= N*(L+M+1)+L*M)
 * @param[out] a State matrix, dimension (LDA,N), column-major
 * @param[in] lda Leading dimension of A (>= max(1,N))
 * @param[out] b Input matrix, dimension (LDB,M), column-major
 * @param[in] ldb Leading dimension of B (>= max(1,N))
 * @param[out] c Output matrix, dimension (LDC,N), column-major
 * @param[in] ldc Leading dimension of C (>= max(1,L))
 * @param[out] d Feedthrough matrix, dimension (LDD,M), column-major
 * @param[in] ldd Leading dimension of D (>= max(1,L))
 * @param[out] x0 Initial state vector, dimension (N)
 * @param[out] dwork Workspace array, dimension (LDWORK)
 * @param[in] ldwork Length of DWORK (>= N*(N+L+1))
 * @param[out] info Exit code: 0 = success, <0 = invalid parameter
 */
void tb01vy(const char* apply, i32 n, i32 m, i32 l, const f64* theta,
            i32 ltheta, f64* a, i32 lda, f64* b, i32 ldb, f64* c, i32 ldc,
            f64* d, i32 ldd, f64* x0, f64* dwork, i32 ldwork, i32* info);

/**
 * @brief Output sequence of linear time-invariant open-loop system.
 *
 * Computes output sequence y(1),...,y(NY) of discrete-time state-space model
 * with system matrix S = [A B; C D] given initial state x(1) and input
 * sequence u(1),...,u(NY).
 *
 * Implements: [x(k+1); y(k)] = S * [x(k); u(k)] for k = 1,...,NY
 *
 * @param[in] n Order of matrix A (n >= 0)
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] p Number of system outputs (p >= 0)
 * @param[in] ny Number of output vectors to compute (ny >= 0)
 * @param[in] s System matrix, dimension (lds,n+m), column-major
 * @param[in] lds Leading dimension of S (lds >= max(1,n+p))
 * @param[in] u Input sequence, dimension (ldu,m), row u(k) contains u(k)'
 * @param[in] ldu Leading dimension of U (ldu >= max(1,ny))
 * @param[in,out] x State vector, dimension (n). On entry: x(1). On exit: x(ny+1)
 * @param[out] y Output sequence, dimension (ldy,p), row y(k) contains y(k)'
 * @param[in] ldy Leading dimension of Y (ldy >= max(1,ny))
 * @param[out] dwork Workspace array, dimension (ldwork)
 * @param[in] ldwork Length of dwork (ldwork >= 2*n+m+p if m>0, n+p if m=0, 0 otherwise)
 * @param[out] info Exit code (0 = success, <0 = invalid parameter)
 */
void tf01mx(const i32 n, const i32 m, const i32 p, const i32 ny,
            const f64* s, const i32 lds, const f64* u, const i32 ldu,
            f64* x, f64* y, const i32 ldy, f64* dwork, const i32 ldwork,
            i32* info);

/**
 * @brief Calculate the output of the Wiener system.
 *
 * @param[in] nsmp Number of training samples.
 * @param[in] m Length of each input sample.
 * @param[in] l Length of each output sample.
 * @param[in] ipar Integer parameters (n, nn).
 * @param[in] lipar Length of ipar.
 * @param[in] x Parameter vector (wb, theta).
 * @param[in] lx Length of x.
 * @param[in] u Input samples.
 * @param[in] ldu Leading dimension of u.
 * @param[out] y Simulated output.
 * @param[in] ldy Leading dimension of y.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01ad(i32 nsmp, i32 m, i32 l, i32 *ipar, i32 lipar, f64 *x, i32 lx, 
            f64 *u, i32 ldu, f64 *y, i32 ldy, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Calculate the Jacobian of the Wiener system.
 *
 * @param[in] cjte 'C' to compute J'*e, 'N' to skip.
 * @param[in] nsmp Number of training samples.
 * @param[in] m Length of each input sample.
 * @param[in] l Length of each output sample.
 * @param[in,out] ipar Integer parameters (n, nn).
 * @param[in] lipar Length of ipar.
 * @param[in,out] x Parameter vector.
 * @param[in] lx Length of x.
 * @param[in] u Input samples.
 * @param[in] ldu Leading dimension of u.
 * @param[in] e Error vector (if cjte='C').
 * @param[out] j Jacobian matrix.
 * @param[in,out] ldj Leading dimension of j.
 * @param[out] jte J'*e product (if cjte='C').
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01bd(const char *cjte, i32 nsmp, i32 m, i32 l, i32 *ipar, i32 lipar, 
            f64 *x, i32 lx, f64 *u, i32 ldu, f64 *e, f64 *j, i32 *ldj, 
            f64 *jte, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Error function for Wiener system identification (FCN for MD03BD).
 *
 * @param[in,out] iflag Integer indicating the action to be performed.
 * @param[in] nsmp Number of training samples.
 * @param[in] n Number of variables.
 * @param[in,out] ipar Integer parameters.
 * @param[in] lipar Length of ipar.
 * @param[in] z Input samples.
 * @param[in] ldz Leading dimension of Z.
 * @param[in] y Output samples.
 * @param[in] ldy Leading dimension of Y.
 * @param[in] x Current estimate of parameters.
 * @param[out] nfevl Number of function evaluations.
 * @param[out] e Error vector.
 * @param[out] j Jacobian matrix.
 * @param[in,out] ldj Leading dimension of J.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01be(i32 *iflag, i32 nsmp, i32 n, i32 *ipar, i32 lipar, 
            f64 *z, i32 ldz, f64 *y, i32 ldy, f64 *x, 
            i32 *nfevl, f64 *e, f64 *j, i32 *ldj, 
            f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Error function for Wiener system identification (Full parameter optimization).
 *
 * @param[in,out] iflag Integer indicating the action to be performed.
 * @param[in] nfun Number of functions.
 * @param[in] lx Number of variables.
 * @param[in,out] ipar Integer parameters.
 * @param[in] lipar Length of ipar.
 * @param[in] u Input samples.
 * @param[in] ldu Leading dimension of U.
 * @param[in] y Output samples.
 * @param[in] ldy Leading dimension of Y.
 * @param[in] x Current estimate of parameters.
 * @param[out] nfevl Number of function evaluations.
 * @param[out] e Error vector.
 * @param[out] j Jacobian matrix.
 * @param[in,out] ldj Leading dimension of J.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01bf(i32 *iflag, i32 nfun, i32 lx, i32 *ipar, i32 lipar, 
            f64 *u, i32 ldu, f64 *y, i32 ldy, f64 *x, 
            i32 *nfevl, f64 *e, f64 *j, i32 *ldj, 
            f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Solve linear system J*x = b, D*x = 0 in least squares sense.
 *
 * @param[in] cond Condition estimation mode.
 * @param[in] n Order of matrix R.
 * @param[in] ipar Integer parameters.
 * @param[in] lipar Length of ipar.
 * @param[in,out] r Matrix R.
 * @param[in] ldr Leading dimension of R.
 * @param[in] ipvt Permutation matrix.
 * @param[in] diag Diagonal scaling.
 * @param[in] qtb Q'*b.
 * @param[in,out] ranks Ranks.
 * @param[out] x Solution.
 * @param[in] tol Tolerance.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01bq(const char *cond, i32 n, const i32 *ipar, i32 lipar, f64 *r, i32 ldr, 
            const i32 *ipvt, const f64 *diag, const f64 *qtb, i32 *ranks, 
            f64 *x, f64 *tol, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Compute Levenberg-Marquardt parameter for Wiener system.
 *
 * @param[in] cond Condition estimation mode.
 * @param[in] n Order of matrix R.
 * @param[in] ipar Integer parameters.
 * @param[in] lipar Length of ipar.
 * @param[in,out] r Matrix R.
 * @param[in] ldr Leading dimension of R.
 * @param[in] ipvt Permutation matrix.
 * @param[in] diag Diagonal scaling.
 * @param[in] qtb Q'*b.
 * @param[in] delta Trust region radius.
 * @param[in,out] par LM parameter.
 * @param[in,out] ranks Ranks.
 * @param[out] x Solution.
 * @param[out] rx Residual.
 * @param[in] tol Tolerance.
 * @param[out] dwork Workspace.
 * @param[in] ldwork Length of dwork.
 * @param[out] info Exit code.
 */
void nf01bp(const char *cond, i32 n, const i32 *ipar, i32 lipar, f64 *r, i32 ldr,
            const i32 *ipvt, const f64 *diag, const f64 *qtb, f64 delta,
            f64 *par, i32 *ranks, f64 *x, f64 *rx, f64 tol, f64 *dwork,
            i32 ldwork, i32 *info);

/**
 * @brief Triangular symmetric rank-k update.
 *
 * Computes either the upper or lower triangular part of:
 *   R = alpha*R + beta*op(A)*B  (SIDE='L')
 *   R = alpha*R + beta*B*op(A)  (SIDE='R')
 *
 * where op(A) = A or A', and only the specified triangle is computed.
 *
 * @param[in] side 'L' for left (R = alpha*R + beta*op(A)*B)
 *                 'R' for right (R = alpha*R + beta*B*op(A))
 * @param[in] uplo 'U' for upper triangle, 'L' for lower triangle
 * @param[in] trans 'N' for op(A)=A, 'T'/'C' for op(A)=A'
 * @param[in] m Order of matrix R
 * @param[in] n Dimension for product:
 *              SIDE='L': rows of B, columns of op(A)
 *              SIDE='R': rows of op(A), columns of B
 * @param[in] alpha Scalar multiplier for R
 * @param[in] beta Scalar multiplier for product
 * @param[in,out] r On entry: m-by-m matrix R (triangle only)
 *                  On exit: updated R (triangle only)
 * @param[in] ldr Leading dimension of R, >= max(1,m)
 * @param[in] a Matrix A with dimensions:
 *              SIDE='L', TRANS='N': m-by-n
 *              SIDE='L', TRANS='T': n-by-m
 *              SIDE='R', TRANS='N': n-by-m
 *              SIDE='R', TRANS='T': m-by-n
 * @param[in] lda Leading dimension of A
 * @param[in] b Matrix B with dimensions:
 *              SIDE='L': n-by-m
 *              SIDE='R': m-by-n
 * @param[in] ldb Leading dimension of B
 * @return 0 on success, -i if argument i had an illegal value
 *
 * @note Main application: computing symmetric updates where B = X*op(A)'
 *       or B = op(A)'*X with X symmetric.
 */
i32 slicot_mb01rx(
    char side,
    char uplo,
    char trans,
    i32 m,
    i32 n,
    f64 alpha,
    f64 beta,
    f64 *r,
    i32 ldr,
    const f64 *a,
    i32 lda,
    const f64 *b,
    i32 ldb
);

/**
 * @brief Apply Householder reflector H to matrix [A; B] from the left.
 *
 * Applies elementary reflector H to (m+1)-by-n matrix C = [A; B],
 * where A has one row:
 *
 *     H = I - tau * u * u',  u = [1; v]
 *
 * Computes C := H * C.
 *
 * Uses inline code for order < 11, BLAS for larger orders.
 *
 * @param[in] m Number of rows of matrix B, m >= 0
 * @param[in] n Number of columns, n >= 0
 * @param[in] v Householder vector v of dimension m
 * @param[in] tau Scalar factor tau (if tau=0, H is identity)
 * @param[in,out] a On entry: 1-by-n matrix A
 *                  On exit: updated first row of H*C
 * @param[in] lda Leading dimension of A, lda >= 1
 * @param[in,out] b On entry: m-by-n matrix B
 *                  On exit: updated last m rows of H*C
 * @param[in] ldb Leading dimension of B, ldb >= max(1,m)
 * @param[out] dwork Workspace of dimension n (not referenced if m+1 < 11)
 *
 * @note Based on LAPACK's DLARFX and DLATZM with special structure optimization.
 */
void SLC_MB04OY(i32 m, i32 n, const f64* v, f64 tau,
                f64* a, i32 lda, f64* b, i32 ldb, f64* dwork);

/**
 * @brief Apply Householder reflector H to matrix [A B] from the right.
 *
 * Applies elementary reflector H to m-by-(n+1) matrix C = [A B],
 * where A has one column:
 *
 *     H = I - tau * u * u',  u = [1; v]
 *
 * Computes C := C * H.
 *
 * Uses inline code for order < 11, BLAS for larger orders.
 *
 * @param[in] m Number of rows of matrices A and B, m >= 0
 * @param[in] n Number of columns of matrix B, n >= 0
 * @param[in] v Householder vector v, dimension (1+(n-1)*abs(incv))
 * @param[in] incv Increment between elements of v, incv != 0
 * @param[in] tau Scalar factor tau (if tau=0, H is identity)
 * @param[in,out] a On entry: m-by-1 matrix A
 *                  On exit: updated first column of C*H
 * @param[in] lda Leading dimension of A, lda >= max(1,m)
 * @param[in,out] b On entry: m-by-n matrix B
 *                  On exit: updated last n columns of C*H
 * @param[in] ldb Leading dimension of B, ldb >= max(1,m)
 * @param[out] dwork Workspace of dimension m (not referenced if n+1 < 11)
 *
 * @note Based on LAPACK's DLARFX and DLATZM with special structure optimization.
 */
void SLC_MB04NY(i32 m, i32 n, const f64* v, i32 incv, f64 tau,
                f64* a, i32 lda, f64* b, i32 ldb, f64* dwork);

/**
 * @brief RQ factorization of special structured block matrix.
 *
 * Calculates an RQ factorization of the first block row and applies
 * the orthogonal transformations (from the right) to the second block row:
 *
 *     [ A   R ]        [ 0   R_new ]
 *     [       ] * Q' = [           ]
 *     [ C   B ]        [ C_new B_new ]
 *
 * where R and R_new are upper triangular. Matrix A can be full or
 * upper trapezoidal/triangular.
 *
 * @param[in] uplo 'U' = A is upper trapezoidal, 'F' = A is full
 * @param[in] n Order of matrices R and R_new, n >= 0
 * @param[in] m Number of rows of matrices B and C, m >= 0
 * @param[in] p Number of columns of matrices A and C, p >= 0
 * @param[in,out] r On entry: n-by-n upper triangular matrix R
 *                  On exit: n-by-n upper triangular matrix R_new
 * @param[in] ldr Leading dimension of R, ldr >= max(1,n)
 * @param[in,out] a On entry: n-by-p matrix A (full or upper trapezoidal)
 *                  On exit: Householder vectors in corresponding positions
 * @param[in] lda Leading dimension of A, lda >= max(1,n)
 * @param[in,out] b On entry: m-by-n matrix B
 *                  On exit: transformed matrix B_new
 * @param[in] ldb Leading dimension of B, ldb >= max(1,m)
 * @param[in,out] c On entry: m-by-p matrix C
 *                  On exit: transformed matrix C_new
 * @param[in] ldc Leading dimension of C, ldc >= max(1,m)
 * @param[out] tau Scalar factors of elementary reflectors, dimension (n)
 * @param[out] dwork Workspace, dimension max(n-1, m)
 */
void SLC_MB04ND(const char* uplo, i32 n, i32 m, i32 p,
                f64* r, i32 ldr, f64* a, i32 lda,
                f64* b, i32 ldb, f64* c, i32 ldc,
                f64* tau, f64* dwork);

/**
 * @brief Apply orthogonal transformations from MB04ID to matrix C.
 *
 * Overwrites real n-by-m matrix C with Q'*C, Q*C, C*Q', or C*Q, where
 * Q is orthogonal matrix defined as product of k elementary reflectors:
 *
 *     Q = H(1) H(2) ... H(k)
 *
 * as returned by MB04ID. Q is order n if SIDE='L', order m if SIDE='R'.
 * Reflectors stored in special format for lower-left zero triangle structure.
 *
 * @param[in] side 'L' = apply Q or Q' from left, 'R' = from right
 * @param[in] trans 'N' = apply Q (no transpose), 'T' = apply Q' (transpose)
 * @param[in] n Number of rows of matrix C, n >= 0
 * @param[in] m Number of columns of matrix C, m >= 0
 * @param[in] k Number of elementary reflectors, constraints:
 *              n >= k >= 0 if SIDE='L'
 *              m >= k >= 0 if SIDE='R'
 * @param[in] p Order of zero triangle (rows of zero trapezoid), p >= 0
 * @param[in,out] a Reflector storage, dimension (lda,k)
 *                  Row i+1:min(n,n-p-1+i) of column i contains H(i)
 *                  Modified but restored on exit
 * @param[in] lda Leading dimension: lda >= max(1,n) if SIDE='L',
 *                                   lda >= max(1,m) if SIDE='R'
 * @param[in] tau Reflector scalar factors, dimension (k)
 * @param[in,out] c On entry: n-by-m matrix C
 *                  On exit: transformed matrix
 * @param[in] ldc Leading dimension, ldc >= max(1,n)
 * @param[out] dwork Workspace, dimension (ldwork)
 *                   dwork[0] returns optimal ldwork on exit
 * @param[in] ldwork Workspace size:
 *                   ldwork >= max(1,m) if SIDE='L'
 *                   ldwork >= max(1,n) if SIDE='R'
 * @param[out] info Exit code: 0=success, <0=invalid parameter
 */
void mb04iy(
    const char* side,
    const char* trans,
    const i32 n,
    const i32 m,
    const i32 k,
    const i32 p,
    f64* a,
    const i32 lda,
    const f64* tau,
    f64* c,
    const i32 ldc,
    f64* dwork,
    const i32 ldwork,
    i32* info
);

/**
 * @brief Compute singular value decomposition of upper triangular matrix.
 *
 * Computes SVD of N-by-N upper triangular matrix: A = Q*S*P', where Q and P
 * are orthogonal matrices and S is diagonal with non-negative singular values
 * in descending order. Uses bidiagonalization followed by QR algorithm.
 *
 * @param[in] jobq 'V' to compute left singular vectors Q, 'N' otherwise
 * @param[in] jobp 'V' to compute right singular vectors P', 'N' otherwise
 * @param[in] n Order of matrix A (n >= 0)
 * @param[in,out] a DOUBLE PRECISION array, dimension (lda,n)
 *                  In: upper triangular matrix A
 *                  Out: if jobp='V', orthogonal matrix P'; otherwise workspace
 * @param[in] lda Leading dimension (lda >= max(1,n))
 * @param[out] q DOUBLE PRECISION array, dimension (ldq,n)
 *               If jobq='V', contains orthogonal matrix Q (left singular vectors)
 * @param[in] ldq Leading dimension (ldq >= max(1,n) if jobq='V', ldq >= 1 otherwise)
 * @param[out] sv DOUBLE PRECISION array, dimension (n)
 *                Singular values in descending order
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   On exit: dwork[0] = optimal ldwork
 *                   If info > 0: dwork[1:n-1] = unconverged superdiagonals
 * @param[in] ldwork Workspace size (ldwork >= max(1,5*n), or -1 for query)
 * @param[out] info Exit code:
 *                  = 0: success
 *                  < 0: if info = -i, i-th argument invalid
 *                  > 0: QR algorithm failed to converge (info = # unconverged superdiagonals)
 */
i32 mb03ud(char jobq, char jobp, i32 n, f64 *a, i32 lda, f64 *q, i32 ldq,
           f64 *sv, f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief User's confirmation of the system order.
 *
 * Non-interactive version for library use. Validates parameters and allows
 * programmatic modification of the estimated system order.
 *
 * In the original Fortran version, this routine provides interactive user
 * confirmation via terminal I/O. For library use, this version validates
 * parameters and ensures N <= NMAX.
 *
 * @param[in] ns Number of singular values (ns > 0)
 * @param[in] nmax Maximum value of system order (0 <= nmax <= ns)
 * @param[in,out] n On entry: estimated system order (0 <= n <= ns)
 *                  On exit: validated order (n <= nmax)
 * @param[in] sv Singular values array, dimension (ns), descending order
 * @param[out] info Exit code:
 *                  = 0: successful exit
 *                  < 0: if info = -i, the i-th argument had an illegal value
 *
 * @note This routine is typically called by IB01OD for system order validation.
 */
i32 SLC_IB01OY(i32 ns, i32 nmax, i32 *n, const f64 *sv, i32 *info);

/**
 * @brief Estimate system order from Hankel singular values.
 *
 * Estimates the system order based on singular values of the relevant part
 * of the triangular factor from QR factorization of concatenated block
 * Hankel matrices.
 *
 * @param[in] ctrl Control mode:
 *                 'C' = call IB01OY for user confirmation
 *                 'N' = no confirmation
 * @param[in] nobr Number of block rows s in Hankel matrices (nobr > 0)
 * @param[in] l Number of system outputs (l > 0)
 * @param[in] sv Singular values array, dimension (l*nobr), descending order
 * @param[out] n Estimated system order
 * @param[in] tol Tolerance for order estimation:
 *                tol >= 0: n = index of last SV >= tol
 *                tol = 0: default tol = nobr * eps * sv[0]
 *                tol < 0: n = index of largest logarithmic gap
 * @param[out] iwarn Warning indicator:
 *                   0 = no warning
 *                   3 = all SVs zero, n = 0
 * @param[out] info Exit code:
 *                  = 0: success
 *                  < 0: if info = -i, parameter i had illegal value
 * @return info value
 */
i32 SLC_IB01OD(char ctrl, i32 nobr, i32 l, const f64 *sv, i32 *n,
               f64 tol, i32 *iwarn, i32 *info);

/**
 * @brief SVD system order via block Hankel.
 *
 * Computes singular value decomposition (SVD) of triangular factor R from
 * QR factorization of concatenated block Hankel matrices to determine system
 * order. Related preliminary calculations for computing system matrices are
 * also performed.
 *
 * @param[in] meth Subspace identification method:
 *                 'M' = MOESP with past inputs/outputs
 *                 'N' = N4SID algorithm
 * @param[in] jobd MOESP BD computation mode (not relevant for N4SID):
 *                 'M' = compute B,D using MOESP approach
 *                 'N' = don't compute B,D using MOESP
 * @param[in] nobr Number of block rows s (nobr > 0)
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] l Number of system outputs (l > 0)
 * @param[in,out] r DOUBLE PRECISION array, dimension (ldr, 2*(m+l)*nobr)
 *                  In: upper triangular factor R from QR factorization
 *                  Out: processed matrix S for subsequent routines
 * @param[in] ldr Leading dimension of R
 *                MOESP/JOBD='M': ldr >= max(2*(m+l)*nobr, 3*m*nobr)
 *                Otherwise: ldr >= 2*(m+l)*nobr
 * @param[out] sv Singular values array, dimension (l*nobr), descending
 * @param[in] tol Tolerance for rank estimation (N4SID only):
 *                tol > 0: lower bound for reciprocal condition number
 *                tol <= 0: use default m*n*eps
 * @param[out] iwork INTEGER array, dimension ((m+l)*nobr)
 *                   Not referenced for METH='M'
 * @param[out] dwork DOUBLE PRECISION array, dimension (ldwork)
 *                   dwork[0] = optimal ldwork
 *                   For N4SID: dwork[1], dwork[2] = reciprocal cond numbers
 * @param[in] ldwork Length of dwork:
 *                   MOESP/JOBD='M': max((2*m-1)*nobr, (m+l)*nobr, 5*l*nobr)
 *                   MOESP/JOBD='N': 5*l*nobr
 *                   N4SID: 5*(m+l)*nobr+1
 * @param[out] iwarn Warning indicator:
 *                   0 = no warning
 *                   4 = U_f rank-deficient (N4SID)
 *                   5 = r_1 rank-deficient (N4SID)
 * @param[out] info Exit code:
 *                  0 = success
 *                  -i = parameter i had illegal value
 *                  2 = SVD did not converge
 * @return info value
 */
i32 SLC_IB01ND(char meth, char jobd, i32 nobr, i32 m, i32 l,
               f64 *r, i32 ldr, f64 *sv, f64 tol,
               i32 *iwork, f64 *dwork, i32 ldwork,
               i32 *iwarn, i32 *info);

/**
 * @brief Estimate initial state and system matrices B, D.
 *
 * Given (A, C) and input/output trajectories, estimates the system matrices
 * B and D, and optionally the initial state x(0), for a discrete-time LTI
 * system: x(k+1) = A*x(k) + B*u(k), y(k) = C*x(k) + D*u(k).
 * Matrix A is assumed to be in real Schur form.
 *
 * @param[in] jobx0 'X' to compute initial state, 'N' if x(0) known to be zero
 * @param[in] job 'B' to compute B only (D known zero), 'D' to compute B and D
 * @param[in] n System order (n >= 0)
 * @param[in] m Number of inputs (m >= 0)
 * @param[in] l Number of outputs (l > 0)
 * @param[in] nsmp Number of samples (rows of U and Y)
 * @param[in] a N-by-N state matrix A in real Schur form, dimension (lda,n)
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in] c L-by-N output matrix C, dimension (ldc,n)
 * @param[in] ldc Leading dimension of C (ldc >= l)
 * @param[in,out] u NSMP-by-M input data, dimension (ldu,m)
 *                  If JOB='D': on exit contains QR factorization details
 * @param[in] ldu Leading dimension of U (ldu >= max(1,nsmp) if m>0, else >= 1)
 * @param[in] y NSMP-by-L output data, dimension (ldy,l)
 * @param[in] ldy Leading dimension of Y (ldy >= max(1,nsmp))
 * @param[out] x0 Estimated initial state, dimension (n)
 * @param[out] b Estimated N-by-M input matrix B, dimension (ldb,m)
 * @param[in] ldb Leading dimension of B (ldb >= n if n>0 and m>0, else >= 1)
 * @param[out] d Estimated L-by-M direct transmission matrix D, dimension (ldd,m)
 * @param[in] ldd Leading dimension of D (ldd >= l if m>0 and JOB='D', else >= 1)
 * @param[in] tol Tolerance for rank estimation (tol <= 0 uses machine precision)
 * @param[out] iwork INTEGER workspace
 * @param[out] dwork DOUBLE PRECISION workspace
 *                   dwork[0] = optimal ldwork
 *                   dwork[1] = rcond of W2 triangular factor
 *                   dwork[2] = rcond of U triangular factor (if JOB='D' and m>0)
 * @param[in] ldwork Workspace size
 * @param[out] iwarn Warning indicator (4 = rank-deficient coefficient matrix)
 * @param[out] info Exit code (0=success, -i=param i invalid, 2=SVD failed)
 */
void slicot_ib01qd(
    const char* jobx0, const char* job,
    i32 n, i32 m, i32 l, i32 nsmp,
    const f64* a, i32 lda,
    const f64* c, i32 ldc,
    f64* u, i32 ldu,
    const f64* y, i32 ldy,
    f64* x0,
    f64* b, i32 ldb,
    f64* d, i32 ldd,
    f64 tol,
    i32* iwork, f64* dwork, i32 ldwork,
    i32* iwarn, i32* info
);

/**
 * @brief System identification driver - MOESP/N4SID preprocessing and order estimation.
 *
 * Preprocesses input-output data for estimating state-space matrices and finds
 * an estimate of the system order using MOESP or N4SID subspace identification.
 * This driver calls IB01MD (R factor), IB01ND (SVD), IB01OD (order estimation).
 *
 * @param[in] meth Method: 'M' = MOESP, 'N' = N4SID
 * @param[in] alg Algorithm: 'C' = Cholesky, 'F' = Fast QR, 'Q' = QR
 * @param[in] jobd MOESP B/D mode: 'M' = compute via MOESP, 'N' = don't (N4SID: not used)
 * @param[in] batch Processing: 'F' = first, 'I' = intermediate, 'L' = last, 'O' = one block
 * @param[in] conct Connection: 'C' = connected, 'N' = not connected (unused if BATCH='O')
 * @param[in] ctrl Confirmation: 'C' = user confirmation via IB01OY, 'N' = no confirmation
 * @param[in] nobr Number of block rows (nobr > 0)
 * @param[in] m Number of inputs (m >= 0)
 * @param[in] l Number of outputs (l > 0)
 * @param[in] nsmp Number of samples (nsmp >= 2*nobr for sequential,
 *                 nsmp >= 2*(m+l+1)*nobr - 1 for non-sequential)
 * @param[in] u NSMP-by-M input data, dimension (ldu,m)
 * @param[in] ldu Leading dimension of U (ldu >= nsmp if m>0, else >= 1)
 * @param[in] y NSMP-by-L output data, dimension (ldy,l)
 * @param[in] ldy Leading dimension of Y (ldy >= nsmp)
 * @param[out] n Estimated system order
 * @param[out] r Upper triangular R/S factor, dimension (ldr, 2*(m+l)*nobr)
 * @param[in] ldr Leading dimension of R (ldr >= max(2*(m+l)*nobr, 3*m*nobr) for MOESP/JOBD='M')
 * @param[out] sv Singular values, dimension (l*nobr)
 * @param[in] rcond Rank tolerance for N4SID (rcond > 0, or <= 0 for default)
 * @param[in] tol Order estimation tolerance:
 *                tol >= 0: n = last SV >= tol; tol = 0: default tol = nobr*eps*sv[0]
 *                tol < 0: n = index of largest logarithmic gap
 * @param[in,out] iwork INTEGER workspace, dimension >= max(3,(m+l)*nobr) for N4SID
 *                      For sequential: iwork[0:2] preserves state
 * @param[out] dwork DOUBLE PRECISION workspace
 * @param[in] ldwork Workspace size (see IB01MD for requirements)
 * @param[out] iwarn Warning: 0=none, 1=100 cycles, 2=fast failed, 3=all SV zero,
 *                   4=U_f rank-deficient, 5=r_1 rank-deficient
 * @param[out] info Exit code: 0=success, -i=param i invalid, 1=fast failed, 2=SVD failed
 */
void ib01ad(const char *meth, const char *alg, const char *jobd,
            const char *batch, const char *conct, const char *ctrl,
            i32 nobr, i32 m, i32 l, i32 nsmp,
            const f64 *u, i32 ldu, const f64 *y, i32 ldy,
            i32 *n, f64 *r, i32 ldr, f64 *sv, f64 rcond, f64 tol,
            i32 *iwork, f64 *dwork, i32 ldwork,
            i32 *iwarn, i32 *info);

/**
 * @brief Upper triangular factor R of concatenated block Hankel matrices.
 *
 * Constructs the upper triangular factor R of the concatenated block
 * Hankel matrices using input-output data. Used in subspace identification
 * methods (MOESP and N4SID). Data can optionally be processed sequentially.
 *
 * For MOESP (METH='M'): H = [Uf' Up' Y']
 * For N4SID (METH='N'): H = [U' Y']
 *
 * @param[in] meth Method: 'M' = MOESP, 'N' = N4SID
 * @param[in] alg Algorithm: 'C' = Cholesky, 'F' = Fast QR, 'Q' = QR
 * @param[in] batch Processing mode: 'F' = first, 'I' = intermediate,
 *                  'L' = last, 'O' = one block only
 * @param[in] conct Connection: 'C' = connected blocks, 'N' = not connected
 * @param[in] nobr Number of block rows s in Hankel matrices (nobr > 0)
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] l Number of system outputs (l > 0)
 * @param[in] nsmp Number of samples (nsmp >= 2*nobr for sequential,
 *                 nsmp >= 2*(m+l+1)*nobr - 1 for non-sequential)
 * @param[in] u NSMP-by-M input data, dimension (ldu,m)
 * @param[in] ldu Leading dimension of U (ldu >= nsmp if m>0, else >= 1)
 * @param[in] y NSMP-by-L output data, dimension (ldy,l)
 * @param[in] ldy Leading dimension of Y (ldy >= nsmp)
 * @param[in,out] r On exit: 2*(m+l)*nobr-by-2*(m+l)*nobr upper triangular R
 *                  On entry for sequential: previous R matrix
 * @param[in] ldr Leading dimension of R (ldr >= 2*(m+l)*nobr)
 * @param[in,out] iwork INTEGER workspace, dimension >= max(3,m+l)
 *                      For sequential: iwork[0:2] preserves state between calls
 * @param[out] dwork DOUBLE PRECISION workspace
 *                   dwork[0] = optimal ldwork on exit
 * @param[in] ldwork Workspace size (use -1 for query)
 * @param[out] iwarn Warning: 0=none, 1=100 cycles exhausted, 2=fast alg failed
 * @param[out] info Exit code: 0=success, -i=param i invalid, 1=fast alg failed
 */
void ib01md(const char *meth, const char *alg, const char *batch,
            const char *conct, i32 nobr, i32 m, i32 l, i32 nsmp,
            const f64 *u, i32 ldu, const f64 *y, i32 ldy,
            f64 *r, i32 ldr, i32 *iwork, f64 *dwork, i32 ldwork,
            i32 *iwarn, i32 *info);

/**
 * @brief Fast QR factorization for block Hankel matrices.
 *
 * Constructs the upper triangular factor R of concatenated block Hankel
 * matrices using input-output data via fast QR based on displacement rank.
 * This is a helper routine called by IB01MD when ALG='F'.
 *
 * @param[in] meth Method: 'M' = MOESP, 'N' = N4SID
 * @param[in] batch Processing mode: 'F' = first, 'I' = intermediate,
 *                  'L' = last, 'O' = one block only
 * @param[in] conct Connection: 'C' = connected blocks, 'N' = not connected
 * @param[in] nobr Number of block rows s in Hankel matrices (nobr > 0)
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] l Number of system outputs (l > 0)
 * @param[in] nsmp Number of samples
 * @param[in] u NSMP-by-M input data, dimension (ldu,m)
 * @param[in] ldu Leading dimension of U
 * @param[in] y NSMP-by-L output data, dimension (ldy,l)
 * @param[in] ldy Leading dimension of Y
 * @param[out] r Upper triangular R factor, dimension (ldr,2*(m+l)*nobr)
 * @param[in] ldr Leading dimension of R (ldr >= 2*(m+l)*nobr)
 * @param[in,out] iwork INTEGER workspace
 * @param[out] dwork DOUBLE PRECISION workspace
 * @param[in] ldwork Workspace size (use -1 for query)
 * @param[out] iwarn Warning indicator
 * @param[out] info Exit code: 0=success, -i=param i invalid, 1=H'H not pos def
 */
void ib01my(const char *meth, const char *batch, const char *conct,
            i32 nobr, i32 m, i32 l, i32 nsmp,
            const f64 *u, i32 ldu, const f64 *y, i32 ldy,
            f64 *r, i32 ldr, i32 *iwork, f64 *dwork, i32 ldwork,
            i32 *iwarn, i32 *info);

/**
 * @brief Estimate system matrices from R factor (subspace identification).
 *
 * Estimates state-space matrices A, C, B, D from the R factor produced by
 * IB01MD. Optionally computes covariance matrices for Kalman gain.
 *
 * @param[in] meth Method: 'M' = MOESP, 'N' = N4SID
 * @param[in] job Matrices to compute: 'A' = all, 'C' = A,C only,
 *                'B' = B only, 'D' = B,D only
 * @param[in] jobcv Covariances: 'C' = compute, 'N' = do not compute
 * @param[in] nobr Block rows (nobr > 1)
 * @param[in] n System order (0 < n < nobr)
 * @param[in] m Number of inputs (m >= 0)
 * @param[in] l Number of outputs (l > 0)
 * @param[in] nsmpl Number of samples (nsmpl >= 2*(m+l)*nobr if jobcv='C')
 * @param[in,out] r R factor from IB01MD, dimension (ldr, 2*(m+l)*nobr)
 * @param[in] ldr Leading dimension of R (ldr >= 2*(m+l)*nobr)
 * @param[in,out] a N-by-N state matrix, dimension (lda,n)
 * @param[in] lda Leading dimension of A
 * @param[in,out] c L-by-N output matrix, dimension (ldc,n)
 * @param[in] ldc Leading dimension of C
 * @param[out] b N-by-M input matrix, dimension (ldb,m)
 * @param[in] ldb Leading dimension of B
 * @param[out] d L-by-M feedthrough matrix, dimension (ldd,m)
 * @param[in] ldd Leading dimension of D
 * @param[out] q N-by-N state covariance, dimension (ldq,n)
 * @param[in] ldq Leading dimension of Q
 * @param[out] ry L-by-L output covariance, dimension (ldry,l)
 * @param[in] ldry Leading dimension of RY
 * @param[out] s N-by-L state-output cross-covariance, dimension (lds,l)
 * @param[in] lds Leading dimension of S
 * @param[out] o L*nobr-by-N extended observability matrix, dimension (ldo,n)
 * @param[in] ldo Leading dimension of O
 * @param[in] tol Tolerance for rank estimation
 * @param[out] iwork INTEGER workspace
 * @param[out] dwork DOUBLE PRECISION workspace
 * @param[in] ldwork Workspace size
 * @param[out] iwarn Warning indicator
 * @param[out] info Exit code
 */
void ib01pd(const char *meth, const char *job, const char *jobcv,
            i32 nobr, i32 n, i32 m, i32 l, i32 nsmpl,
            f64 *r, i32 ldr, f64 *a, i32 lda, f64 *c, i32 ldc,
            f64 *b, i32 ldb, f64 *d, i32 ldd, f64 *q, i32 ldq,
            f64 *ry, i32 ldry, f64 *s, i32 lds, f64 *o, i32 ldo,
            f64 tol, i32 *iwork, f64 *dwork, i32 ldwork,
            i32 *iwarn, i32 *info);

/**
 * @brief Estimate state-space matrices from N4SID/MOESP triangular factor.
 *
 * Estimates system matrices (A,C,B,D), optionally noise covariance matrices
 * (Q,Ry,S), and Kalman gain K from the triangular factor R computed by IB01AD.
 * Supports N4SID, MOESP, or combined subspace identification methods.
 *
 * @param[in] meth Method: 'M' MOESP, 'N' N4SID, 'C' combined
 * @param[in] job Job: 'A' all matrices, 'C' A,C only, 'B' B,D only, 'D' D only
 * @param[in] jobck Covariance: 'K' Kalman gain, 'C' covariances, 'N' neither
 * @param[in] nobr Number of block rows s >= 2
 * @param[in] n System order (0 < n < nobr)
 * @param[in] m Number of inputs (m >= 0)
 * @param[in] l Number of outputs (l > 0)
 * @param[in] nsmpl Number of samples (nsmpl >= 2*(m+l)*nobr for covariances)
 * @param[in,out] r Triangular factor from IB01AD, dimension (ldr, 2*(m+l)*nobr)
 * @param[in] ldr Leading dimension of R
 * @param[out] a N-by-N state matrix, dimension (lda,n)
 * @param[in] lda Leading dimension of A
 * @param[out] c L-by-N output matrix, dimension (ldc,n)
 * @param[in] ldc Leading dimension of C
 * @param[out] b N-by-M input matrix, dimension (ldb,m)
 * @param[in] ldb Leading dimension of B
 * @param[out] d L-by-M feedthrough matrix, dimension (ldd,m)
 * @param[in] ldd Leading dimension of D
 * @param[out] q N-by-N state covariance, dimension (ldq,n)
 * @param[in] ldq Leading dimension of Q
 * @param[out] ry L-by-L output covariance, dimension (ldry,l)
 * @param[in] ldry Leading dimension of RY
 * @param[out] s N-by-L state-output cross-covariance, dimension (lds,l)
 * @param[in] lds Leading dimension of S
 * @param[out] k N-by-L Kalman gain, dimension (ldk,l)
 * @param[in] ldk Leading dimension of K
 * @param[in] tol Tolerance for rank estimation
 * @param[out] iwork INTEGER workspace
 * @param[out] dwork DOUBLE PRECISION workspace
 * @param[in] ldwork Workspace size
 * @param[out] iwarn Warning indicator
 * @param[out] info Exit code
 */
void ib01bd(const char *meth, const char *job, const char *jobck,
            i32 nobr, i32 n, i32 m, i32 l, i32 nsmpl,
            f64 *r, i32 ldr, f64 *a, i32 lda, f64 *c, i32 ldc,
            f64 *b, i32 ldb, f64 *d, i32 ldd, f64 *q, i32 ldq,
            f64 *ry, i32 ldry, f64 *s, i32 lds, f64 *k, i32 ldk,
            f64 tol, i32 *iwork, f64 *dwork, i32 ldwork, i32 *bwork,
            i32 *iwarn, i32 *info);

void ib01px(const char *job, i32 nobr, i32 n, i32 m, i32 l,
            f64 *uf, i32 lduf, const f64 *un, i32 ldun,
            f64 *ul, i32 ldul, const f64 *pgal, i32 ldpgal,
            const f64 *k, i32 ldk, f64 *r, i32 ldr, f64 *x,
            f64 *b, i32 ldb, f64 *d, i32 ldd, f64 tol,
            i32 *iwork, f64 *dwork, i32 ldwork, i32 *iwarn, i32 *info);

void ib01py(const char *meth, const char *job, i32 nobr, i32 n, i32 m, i32 l,
            i32 rankr1, f64 *ul, i32 ldul, const f64 *r1, i32 ldr1,
            const f64 *tau1, const f64 *pgal, i32 ldpgal,
            f64 *k, i32 ldk, f64 *r, i32 ldr, f64 *h, i32 ldh,
            f64 *b, i32 ldb, f64 *d, i32 ldd, f64 tol,
            i32 *iwork, f64 *dwork, i32 ldwork, i32 *iwarn, i32 *info);

/**
 * @brief Wiener system identification using Levenberg-Marquardt algorithm.
 *
 * Computes parameters for approximating a Wiener system consisting of a
 * linear state-space part and a static nonlinearity (neural network):
 *   x(t+1) = A*x(t) + B*u(t)       (linear state-space)
 *   z(t)   = C*x(t) + D*u(t)
 *   y(t)   = f(z(t), wb(1:L))      (nonlinear function)
 *
 * The parameter vector X = (wb(1),...,wb(L), theta) where wb(i) are neural
 * network weights and theta are linear part parameters in output normal form.
 *
 * @param[in] init Initialization mode:
 *                 'L' = initialize linear part only
 *                 'S' = initialize static nonlinearity only
 *                 'B' = initialize both parts
 *                 'N' = no initialization (use given X)
 * @param[in] nobr Block rows for MOESP/N4SID (used if INIT='L' or 'B')
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] l Number of system outputs (l >= 0, l > 0 if INIT='L' or 'B')
 * @param[in] nsmp Number of input/output samples
 * @param[in,out] n System order. If n < 0 and INIT='L' or 'B', order is estimated.
 * @param[in] nn Number of neurons for nonlinear approximation (nn >= 0)
 * @param[in] itmax1 Max iterations for nonlinear initialization (ignored if INIT='L' or 'N')
 * @param[in] itmax2 Max iterations for whole optimization (itmax2 >= 0)
 * @param[in] nprint Print control (> 0 enables iteration printing)
 * @param[in] u Input samples, dimension (ldu, m)
 * @param[in] ldu Leading dimension of U (ldu >= max(1, nsmp))
 * @param[in] y Output samples, dimension (ldy, l)
 * @param[in] ldy Leading dimension of Y (ldy >= max(1, nsmp))
 * @param[in,out] x Parameter vector, dimension (lx).
 *                  On entry: initial parameters (depending on INIT mode)
 *                  On exit: optimized parameters
 * @param[in,out] lx Length of X. On exit, may be updated if N was auto-detected.
 * @param[in] tol1 Tolerance for nonlinear initialization (tol1 < 0 uses sqrt(eps))
 * @param[in] tol2 Tolerance for whole optimization (tol2 < 0 uses sqrt(eps))
 * @param[out] iwork INTEGER workspace
 *                   On exit: iwork[0]=fcn evals, iwork[1]=jac evals,
 *                   iwork[2]=number of condition numbers in dwork
 * @param[in,out] dwork DOUBLE PRECISION workspace
 *                      On entry: dwork[0:3] = seed for random initialization
 *                      On exit: dwork[0]=opt workspace, dwork[1]=residual,
 *                      dwork[2]=iterations, dwork[3]=final Levenberg factor
 * @param[in] ldwork Workspace size
 * @param[out] iwarn Warning indicator (k*100 + j*10 + i)
 * @param[out] info Exit code (0=success, <0=invalid param, >0=algorithm error)
 */
void ib03bd(
    const char* init,
    i32 nobr, i32 m, i32 l, i32 nsmp,
    i32* n,
    i32 nn, i32 itmax1, i32 itmax2, i32 nprint,
    const f64* u, i32 ldu,
    const f64* y, i32 ldy,
    f64* x, i32* lx,
    f64 tol1, f64 tol2,
    i32* iwork, f64* dwork, i32 ldwork,
    i32* iwarn, i32* info);

void mb01vd(const char *trana, const char *tranb, i32 ma, i32 na, i32 mb, i32 nb,
            f64 alpha, f64 beta, const f64 *a, i32 lda, const f64 *b, i32 ldb,
            f64 *c, i32 ldc, i32 *mc, i32 *nc, i32 *info);

void mb02qy(i32 m, i32 n, i32 nrhs, i32 rank, f64 *a, i32 lda,
            const i32 *jpvt, f64 *b, i32 ldb, f64 *tau,
            f64 *dwork, i32 ldwork, i32 *info);

/**
 * @brief Estimate initial state for discrete-time LTI system.
 *
 * Given system matrices (A,B,C,D) and input/output trajectories, estimates
 * the initial state x(0) of the discrete-time LTI system:
 *   x(k+1) = A*x(k) + B*u(k)
 *   y(k)   = C*x(k) + D*u(k)
 *
 * Matrix A is assumed to be in real Schur form.
 *
 * @param[in] job 'Z' if D matrix is zero, 'N' if D is not zero
 * @param[in] n System order (n >= 0)
 * @param[in] m Number of inputs (m >= 0)
 * @param[in] l Number of outputs (l > 0)
 * @param[in] nsmp Number of samples (nsmp >= n)
 * @param[in] a N-by-N state matrix A in real Schur form, dimension (lda,n)
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in] b N-by-M input matrix B, dimension (ldb,m)
 * @param[in] ldb Leading dimension of B (ldb >= n if n>0 and m>0, else >= 1)
 * @param[in] c L-by-N output matrix C, dimension (ldc,n)
 * @param[in] ldc Leading dimension of C (ldc >= l)
 * @param[in] d L-by-M direct transmission matrix D, dimension (ldd,m)
 * @param[in] ldd Leading dimension of D (ldd >= l if m>0 and job='N', else >= 1)
 * @param[in] u NSMP-by-M input data U, dimension (ldu,m)
 * @param[in] ldu Leading dimension of U (ldu >= max(1,nsmp) if m>0, else >= 1)
 * @param[in] y NSMP-by-L output data Y, dimension (ldy,l)
 * @param[in] ldy Leading dimension of Y (ldy >= max(1,nsmp))
 * @param[out] x0 Estimated initial state, dimension (n)
 * @param[in] tol Tolerance for rank estimation (tol <= 0 uses machine precision)
 * @param[out] iwork INTEGER workspace, dimension (n)
 * @param[out] dwork DOUBLE PRECISION workspace
 *                   dwork[0] = optimal ldwork
 *                   dwork[1] = reciprocal condition number of triangular factor
 * @param[in] ldwork Workspace size
 * @param[out] iwarn Warning indicator (4 = rank-deficient coefficient matrix)
 * @param[out] info Exit code (0=success, -i=param i invalid, 2=SVD failed)
 */
void slicot_ib01rd(
    const char* job,
    i32 n, i32 m, i32 l, i32 nsmp,
    const f64* a, i32 lda,
    const f64* b, i32 ldb,
    const f64* c, i32 ldc,
    const f64* d, i32 ldd,
    const f64* u, i32 ldu,
    const f64* y, i32 ldy,
    f64* x0,
    f64 tol,
    i32* iwork, f64* dwork, i32 ldwork,
    i32* iwarn, i32* info
);

/**
 * @brief Estimate initial state and system matrices B, D (driver routine).
 *
 * Estimates initial state x(0) and optionally B and D for discrete-time LTI:
 *   x(k+1) = A*x(k) + B*u(k)
 *   y(k)   = C*x(k) + D*u(k)
 *
 * Driver routine that:
 * 1. Transforms A to real Schur form via TB01WD (A = V*At*V')
 * 2. Calls IB01QD (COMUSE='C') or IB01RD for estimation
 * 3. Back-transforms results to original coordinates
 *
 * @param[in] jobx0 Initial state computation:
 *                  'X' = compute initial state x(0)
 *                  'N' = do not compute x(0), set to zero
 * @param[in] comuse How to handle B and D matrices:
 *                   'C' = compute B (and D if JOB='D')
 *                   'U' = use given B (and D if JOB='D')
 *                   'N' = do not compute/use B and D
 * @param[in] job Matrix computation extent:
 *                'B' = compute B only (D is zero)
 *                'D' = compute B and D
 * @param[in] n System order (n >= 0)
 * @param[in] m Number of inputs (m >= 0)
 * @param[in] l Number of outputs (l > 0)
 * @param[in] nsmp Number of samples
 * @param[in] a N-by-N state matrix A, dimension (lda,n)
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in,out] b N-by-M input matrix B, dimension (ldb,m)
 *                  If COMUSE='U': on entry, given B matrix
 *                  If COMUSE='C': on exit, estimated B matrix
 * @param[in] ldb Leading dimension of B (ldb >= n if m>0, else >= 1)
 * @param[in] c L-by-N output matrix C, dimension (ldc,n)
 * @param[in] ldc Leading dimension of C (ldc >= l)
 * @param[in,out] d L-by-M feedthrough matrix D, dimension (ldd,m)
 *                  If COMUSE='U' and JOB='D': on entry, given D matrix
 *                  If COMUSE='C' and JOB='D': on exit, estimated D matrix
 * @param[in] ldd Leading dimension of D (ldd >= l if m>0 and JOB='D', else >= 1)
 * @param[in] u NSMP-by-M input data U, dimension (ldu,m)
 * @param[in] ldu Leading dimension of U (ldu >= nsmp if m>0, else >= 1)
 * @param[in] y NSMP-by-L output data Y, dimension (ldy,l)
 * @param[in] ldy Leading dimension of Y (ldy >= nsmp)
 * @param[out] x0 Estimated initial state, dimension (n)
 * @param[out] v N-by-N orthogonal transformation matrix, dimension (ldv,n)
 *               Satisfies A = V * At * V', where At is in Schur form
 * @param[in] ldv Leading dimension of V (ldv >= max(1,n))
 * @param[in] tol Tolerance for rank estimation (tol <= 0 uses default)
 * @param[out] iwork INTEGER workspace
 * @param[out] dwork DOUBLE PRECISION workspace
 *                   dwork[0] = optimal workspace
 *                   dwork[1] = reciprocal condition number
 *                   dwork[2] = rcond of U (if COMUSE='C', JOB='D', M>0)
 * @param[in] ldwork Workspace size
 * @param[out] iwarn Warning: 0=none, 4=rank-deficient, 6=A not stable
 * @param[out] info Exit code: 0=success, 1=Schur failed, -i=param i invalid
 */
void slicot_ib01cd(
    const char* jobx0, const char* comuse, const char* job,
    i32 n, i32 m, i32 l, i32 nsmp,
    const f64* a, i32 lda,
    f64* b, i32 ldb,
    const f64* c, i32 ldc,
    f64* d, i32 ldd,
    f64* u, i32 ldu,
    const f64* y, i32 ldy,
    f64* x0,
    f64* v, i32 ldv,
    f64 tol,
    i32* iwork, f64* dwork, i32 ldwork,
    i32* iwarn, i32* info
);

/**
 * @brief Block symmetric rank-k update (BLAS 3 version of MB01RX).
 *
 * Computes triangular part of matrix formula:
 *   R = alpha*R + beta*op(A)*B  (SIDE='L'), or
 *   R = alpha*R + beta*B*op(A)  (SIDE='R')
 *
 * where alpha, beta are scalars, R is m-by-m, and op(A) is A or A'.
 *
 * @param[in] side 'L' for R=alpha*R+beta*op(A)*B, 'R' for R=alpha*R+beta*B*op(A)
 * @param[in] uplo 'U' for upper triangle, 'L' for lower triangle
 * @param[in] trans 'N' for op(A)=A, 'T'/'C' for op(A)=A'
 * @param[in] m Order of R (m >= 0)
 * @param[in] n Inner dimension (n >= 0)
 * @param[in] alpha Scalar alpha
 * @param[in] beta Scalar beta
 * @param[in,out] r m-by-m matrix R, dimension (ldr,m)
 * @param[in] ldr Leading dimension of R (ldr >= max(1,m))
 * @param[in] a Matrix A
 * @param[in] lda Leading dimension of A
 * @param[in] b Matrix B
 * @param[in] ldb Leading dimension of B
 * @param[out] info Exit code (0=success, <0=invalid parameter)
 *
 * @note Main application: symmetric updates where B = X*op(A)' or B = op(A)'*X
 */
void mb01rb(const char* side, const char* uplo, const char* trans,
            const i32 m, const i32 n, const f64 alpha, const f64 beta,
            f64* r, const i32 ldr, const f64* a, const i32 lda,
            const f64* b, const i32 ldb, i32* info);

/**
 * @brief Scale rows or columns of a matrix by a diagonal matrix.
 *
 * Computes one of:
 *   A := diag(R) * A        (jobs='R', row scaling)
 *   A := A * diag(C)        (jobs='C', column scaling)
 *   A := diag(R) * A * diag(C)  (jobs='B', both)
 *
 * @param[in] jobs Scaling operation: 'R'=row, 'C'=column, 'B'=both
 * @param[in] m Number of rows of A (m >= 0)
 * @param[in] n Number of columns of A (n >= 0)
 * @param[in,out] a M-by-N matrix, dimension (lda,n), scaled on exit
 * @param[in] lda Leading dimension of A (lda >= max(1,m))
 * @param[in] r Row scale factors, dimension (m), not used if jobs='C'
 * @param[in] c Column scale factors, dimension (n), not used if jobs='R'
 */
void mb01sd(const char jobs, const i32 m, const i32 n,
            f64* a, const i32 lda, const f64* r, const f64* c);

/**
 * @brief Product of upper quasi-triangular matrices B := A * B.
 *
 * Computes the matrix product A * B, where A and B are upper quasi-triangular
 * matrices (block upper triangular with 1-by-1 or 2-by-2 diagonal blocks)
 * with the same structure. The result is returned in array B.
 *
 * @param[in] n Order of matrices A and B (n >= 0)
 * @param[in] a N-by-N upper quasi-triangular matrix, dimension (lda,n)
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in,out] b N-by-N upper quasi-triangular matrix with same structure as A.
 *                  On exit: contains the product A * B.
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n))
 * @param[out] dwork Workspace, dimension (n-1)
 * @param[out] info Exit code: 0=success, <0=invalid parameter,
 *                  1=A and B have different structures or are not quasi-triangular
 *
 * @note Useful for computing powers of real Schur form matrices.
 */
void mb01td(const i32 n, const f64* a, const i32 lda,
            f64* b, const i32 ldb, f64* dwork, i32* info);

/**
 * @brief QR factorization of structured block matrix.
 *
 * Calculates QR factorization of first block column and applies orthogonal
 * transformations to second block column:
 *         [ R   B ]   [ R_   B_ ]
 *    Q' * [       ] = [         ]
 *         [ A   C ]   [ 0    C_ ]
 *
 * where R and R_ are upper triangular, A can be full or upper trapezoidal.
 *
 * @param[in] uplo 'U' for A upper trapezoidal/triangular, 'F' for A full
 * @param[in] n Order of R (n >= 0)
 * @param[in] m Number of columns in B, C (m >= 0)
 * @param[in] p Number of rows in A, C (p >= 0)
 * @param[in,out] r n-by-n upper triangular matrix R, dimension (ldr,n)
 * @param[in] ldr Leading dimension of R (ldr >= max(1,n))
 * @param[in,out] a p-by-n matrix A, dimension (lda,n)
 *                  On exit: Householder vectors
 * @param[in] lda Leading dimension of A (lda >= max(1,p))
 * @param[in,out] b n-by-m matrix B, dimension (ldb,m)
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n))
 * @param[in,out] c p-by-m matrix C, dimension (ldc,m)
 * @param[in] ldc Leading dimension of C (ldc >= max(1,p))
 * @param[out] tau Householder scalars, dimension (n)
 * @param[out] dwork Workspace, dimension (max(n-1,m))
 *
 * @note Algorithm is backward stable
 */
void mb04od(const char* uplo, const i32 n, const i32 m, const i32 p,
            f64* r, const i32 ldr, f64* a, const i32 lda,
            f64* b, const i32 ldb, f64* c, const i32 ldc,
            f64* tau, f64* dwork);

/**
 * @brief Apply Householder transformation to block matrix rows.
 *
 * Inline code for applying Householder transformation H = I - tau*u*u' where
 * u = [1; v] to block matrix rows, exploiting structure in MB04OD.
 *
 * @param[in] m Length of vector v (m >= 0)
 * @param[in] n Number of columns to update (n >= 0)
 * @param[in] v Householder vector v, dimension (m)
 * @param[in] tau Householder scalar
 * @param[in,out] c1 First row to update, dimension (n)
 * @param[in] ldc1 Leading dimension of C1
 * @param[in,out] c2 Remaining m rows to update, dimension (ldc2,n)
 * @param[in] ldc2 Leading dimension of C2
 * @param[out] dwork Workspace, dimension (n)
 */
void mb04oy(const i32* m, const i32* n, const f64* v, const f64* tau,
            f64* c1, const i32* ldc1, f64* c2, const i32* ldc2, f64* dwork);

/**
 * @brief Riccati preprocessing - convert coupling weight problems to standard form.
 *
 * Computes:
 *   G = B*R^(-1)*B'
 *   A_bar = A - B*R^(-1)*L'
 *   Q_bar = Q - L*R^(-1)*L'
 *
 * where A, B, Q, R, L, and G are N-by-N, N-by-M, N-by-N, M-by-M,
 * N-by-M, and N-by-N matrices, respectively, with Q, R and G symmetric.
 *
 * @param[in] jobg 'G' to compute G, 'N' to not compute G
 * @param[in] jobl 'Z' if L is zero, 'N' if L is nonzero
 * @param[in] fact 'N' R unfactored, 'C' Cholesky factor, 'U' UdU'/LdL' factor
 * @param[in] uplo 'U' upper triangle stored, 'L' lower triangle stored
 * @param[in] n Order of A, Q, G; rows of B, L (n >= 0)
 * @param[in] m Order of R; columns of B, L (m >= 0)
 * @param[in,out] a N-by-N matrix A (if jobl='N'), dimension (lda,n)
 *                  On exit: A_bar = A - B*R^(-1)*L'
 * @param[in] lda Leading dimension of A (lda >= max(1,n) if jobl='N', else >= 1)
 * @param[in,out] b N-by-M matrix B, dimension (ldb,m)
 *                  On exit if oufact=1: B*chol(R)^(-1)
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n))
 * @param[in,out] q N-by-N symmetric matrix Q (if jobl='N'), dimension (ldq,n)
 *                  On exit: Q_bar = Q - L*R^(-1)*L'
 * @param[in] ldq Leading dimension of Q (ldq >= max(1,n) if jobl='N', else >= 1)
 * @param[in,out] r M-by-M symmetric matrix R, dimension (ldr,m)
 *                  On exit if oufact=1: Cholesky factor
 *                  On exit if oufact=2: UdU'/LdL' factors
 * @param[in] ldr Leading dimension of R (ldr >= max(1,m))
 * @param[in,out] l N-by-M matrix L (if jobl='N'), dimension (ldl,m)
 *                  On exit if oufact=1: L*chol(R)^(-1)
 * @param[in] ldl Leading dimension of L (ldl >= max(1,n) if jobl='N', else >= 1)
 * @param[in,out] ipiv Pivot indices for UdU'/LdL' (dimension m)
 * @param[out] oufact 0=no factorization (m=0), 1=Cholesky, 2=UdU'/LdL'
 * @param[out] g N-by-N matrix G = B*R^(-1)*B' (if jobg='G'), dimension (ldg,n)
 * @param[in] ldg Leading dimension of G (ldg >= max(1,n) if jobg='G', else >= 1)
 * @param[out] iwork Integer workspace, dimension (m). Not referenced if fact='C'/'U'.
 * @param[out] dwork Double workspace, dimension (ldwork)
 *                   On exit: dwork[0]=optimal ldwork, dwork[1]=rcond (if fact='N')
 * @param[in] ldwork Workspace size. Required sizes depend on fact, jobg, jobl.
 * @param[out] info 0=success, <0=invalid param, 1..m=singular d factor, m+1=R singular
 */
void sb02mt(
    const char* jobg,
    const char* jobl,
    const char* fact,
    const char* uplo,
    const i32 n,
    const i32 m,
    f64* a,
    const i32 lda,
    f64* b,
    const i32 ldb,
    f64* q,
    const i32 ldq,
    f64* r,
    const i32 ldr,
    f64* l,
    const i32 ldl,
    i32* ipiv,
    i32* oufact,
    f64* g,
    const i32 ldg,
    i32* iwork,
    f64* dwork,
    const i32 ldwork,
    i32* info);

/**
 * @brief Optimal state feedback matrix for optimal control problem.
 *
 * Computes the optimal feedback matrix F for the problem of optimal control:
 *
 *   F = (R + B'XB)^(-1) (B'XA + L')   [discrete-time, DICO='D']
 *   F = R^(-1) (B'X + L')             [continuous-time, DICO='C']
 *
 * where A is N-by-N, B is N-by-M, L is N-by-M, R and X are M-by-M and N-by-N
 * symmetric matrices respectively.
 *
 * @param[in] dico 'D' for discrete-time, 'C' for continuous-time
 * @param[in] fact Specifies how R is given:
 *                 'N' = R is unfactored
 *                 'D' = R contains P-by-M matrix D where R = D'D
 *                 'C' = R contains Cholesky factor
 *                 'U' = R contains UdU'/LdL' factorization (continuous only)
 * @param[in] uplo 'U' = upper triangle stored, 'L' = lower triangle
 * @param[in] jobl 'Z' = L is zero, 'N' = L is nonzero
 * @param[in] n Order of matrices A and X (n >= 0)
 * @param[in] m Number of system inputs (m >= 0)
 * @param[in] p Number of rows of D (fact='D' only, p >= m for continuous)
 * @param[in] a N-by-N state matrix A (discrete only), dimension (lda,n)
 * @param[in] lda Leading dimension of A (lda >= max(1,n) if discrete, else >= 1)
 * @param[in,out] b N-by-M input matrix B, dimension (ldb,m)
 *                  May be modified on exit for discrete-time with factored R
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n))
 * @param[in,out] r M-by-M symmetric input weighting matrix R, dimension (ldr,m)
 *                  On exit: Cholesky or UdU'/LdL' factorization
 *                  For fact='D': P-by-M matrix D on entry
 * @param[in] ldr Leading dimension of R (ldr >= max(1,m), or max(1,m,p) for fact='D')
 * @param[in,out] ipiv Pivot indices for UdU'/LdL' factorization, dimension (m)
 *                     Input for fact='U', output for oufact[0]=2
 * @param[in] l N-by-M cross weighting matrix L (if jobl='N'), dimension (ldl,m)
 * @param[in] ldl Leading dimension of L (ldl >= max(1,n) if jobl='N', else >= 1)
 * @param[in,out] x N-by-N Riccati solution matrix X, dimension (ldx,n)
 *                  May be modified for discrete-time with factored R
 * @param[in] ldx Leading dimension of X (ldx >= max(1,n))
 * @param[in] rnorm 1-norm of original R (required for fact='U' only)
 * @param[out] f M-by-N optimal feedback matrix F, dimension (ldf,n)
 * @param[in] ldf Leading dimension of F (ldf >= max(1,m))
 * @param[out] oufact Array of dimension 2:
 *                    oufact[0]: 1=Cholesky of R/R+B'XB, 2=UdU'/LdL'
 *                    oufact[1]: 1=Cholesky of X, 2=spectral (discrete+factored only)
 * @param[out] dwork Workspace, dimension (ldwork)
 *                   On exit: dwork[0]=optimal ldwork, dwork[1]=rcond
 *                   If oufact[1]=2: dwork[2..n+1] contain eigenvalues of X
 * @param[in] ldwork Workspace size:
 *                   fact='U': >= max(2,2*m)
 *                   fact!='U', dico='C': >= max(2,3*m)
 *                   fact='N', dico='D': >= max(2,3*m,n)
 *                   fact!='N', dico='D': >= max(n+3*m+2,4*n+1)
 * @param[out] info 0=success, <0=invalid param, 1..m=singular d factor,
 *                  m+1=R singular, m+2=eigenvalue convergence, m+3=X indefinite
 */
void sb02nd(
    const char* dico,
    const char* fact,
    const char* uplo,
    const char* jobl,
    const i32 n,
    const i32 m,
    const i32 p,
    f64* a,
    const i32 lda,
    f64* b,
    const i32 ldb,
    f64* r,
    const i32 ldr,
    i32* ipiv,
    f64* l,
    const i32 ldl,
    f64* x,
    const i32 ldx,
    const f64 rnorm,
    f64* f,
    const i32 ldf,
    i32* oufact,
    f64* dwork,
    const i32 ldwork,
    i32* info);

/**
 * @brief Solve algebraic Riccati equation using Schur vector method.
 *
 * Solves continuous-time algebraic Riccati equation:
 *     Q + op(A)'*X + X*op(A) - X*G*X = 0                   (DICO='C')
 *
 * or discrete-time algebraic Riccati equation:
 *     Q + op(A)'*X*(I_n + G*X)^(-1)*op(A) - X = 0         (DICO='D')
 *
 * where op(M) = M or M' (transpose), A, G, Q are N-by-N matrices,
 * G and Q are symmetric, and X is the symmetric solution.
 *
 * The matrix G = op(B)*R^(-1)*op(B)' must be provided instead of B and R.
 * Use SB02MT to compute G from B and R.
 *
 * @param[in] job Computation to perform:
 *                'X' = compute solution only
 *                'C' = compute reciprocal condition number only
 *                'E' = compute error bound only
 *                'A' = compute all (solution, condition, error bound)
 * @param[in] dico Problem type:
 *                 'C' = continuous-time
 *                 'D' = discrete-time
 * @param[in] hinv For discrete-time (DICO='D') with JOB='X'/'A':
 *                 'D' = construct symplectic matrix H
 *                 'I' = construct inverse of H
 * @param[in] trana Form of op(A):
 *                  'N' = op(A) = A
 *                  'T'/'C' = op(A) = A'
 * @param[in] uplo Triangle stored for G and Q:
 *                 'U' = upper triangle
 *                 'L' = lower triangle
 * @param[in] scal Scaling strategy (for JOB='X'/'A'):
 *                 'G' = general scaling
 *                 'N' = no scaling
 * @param[in] sort Eigenvalue ordering (for JOB='X'/'A'):
 *                 'S' = stable eigenvalues first
 *                 'U' = unstable eigenvalues first
 * @param[in] fact Schur factorization (for JOB!='X'):
 *                 'F' = T,V contain Schur factors
 *                 'N' = compute Schur factors
 * @param[in] lyapun Lyapunov equation form (for JOB!='X'):
 *                   'O' = original equations
 *                   'R' = reduced equations
 * @param[in] n Order of matrices A, Q, G, X (n >= 0)
 * @param[in] a N-by-N matrix A, dimension (lda,n)
 * @param[in] lda Leading dimension of A
 * @param[in,out] t N-by-N Schur form matrix (for JOB!='X'), dimension (ldt,n)
 * @param[in] ldt Leading dimension of T
 * @param[in,out] v N-by-N orthogonal Schur matrix (for JOB!='X'), dimension (ldv,n)
 * @param[in] ldv Leading dimension of V
 * @param[in,out] g N-by-N symmetric matrix G, dimension (ldg,n)
 * @param[in] ldg Leading dimension of G
 * @param[in,out] q N-by-N symmetric matrix Q, dimension (ldq,n)
 * @param[in] ldq Leading dimension of Q
 * @param[out] x N-by-N symmetric solution X, dimension (ldx,n)
 * @param[in] ldx Leading dimension of X
 * @param[out] sep Separation or scaling factor
 * @param[out] rcond Reciprocal condition number (for JOB='C'/'A')
 * @param[out] ferr Forward error bound (for JOB='E'/'A')
 * @param[out] wr Real parts of eigenvalues, dimension (2*n)
 * @param[out] wi Imaginary parts of eigenvalues, dimension (2*n)
 * @param[out] s 2N-by-2N ordered Schur form, dimension (lds,2*n)
 * @param[in] lds Leading dimension of S
 * @param[out] iwork Integer workspace
 * @param[out] dwork Double workspace
 * @param[in] ldwork Workspace size (>= 5+max(1,4*n*n+8*n) for JOB='X'/'A')
 * @param[out] bwork Logical workspace, dimension (2*n)
 * @param[out] info 0=success, 1=A singular, 2=Schur failed, 3=ordering failed,
 *                  4=not enough stable eigenvalues, 5=linear system singular,
 *                  6=Ac Schur failed, 7=near-equal eigenvalues (warning)
 */
void sb02rd(
    const char* job,
    const char* dico,
    const char* hinv,
    const char* trana,
    const char* uplo,
    const char* scal,
    const char* sort,
    const char* fact,
    const char* lyapun,
    const i32 n,
    f64* a,
    const i32 lda,
    f64* t,
    const i32 ldt,
    f64* v,
    const i32 ldv,
    f64* g,
    const i32 ldg,
    f64* q,
    const i32 ldq,
    f64* x,
    const i32 ldx,
    f64* sep,
    f64* rcond,
    f64* ferr,
    f64* wr,
    f64* wi,
    f64* s,
    const i32 lds,
    i32* iwork,
    f64* dwork,
    const i32 ldwork,
    i32* bwork,
    i32* info);

/**
 * @brief Construct Hamiltonian or symplectic matrix for Riccati equation.
 *
 * For continuous-time (DICO='C'):
 *         ( op(A)   -G    )
 *     S = (               )
 *         (  -Q   -op(A)' )
 *
 * For discrete-time (DICO='D'):
 *                 -1              -1
 *         (  op(A)           op(A)  *G       )
 *     S = (        -1                   -1   )  if HINV='D'
 *         ( Q*op(A)     op(A)' + Q*op(A)  *G )
 *
 * @param[in] dico Problem type: 'C'=continuous, 'D'=discrete
 * @param[in] hinv For DICO='D': 'D'=direct, 'I'=inverse
 * @param[in] trana Form of op(A): 'N'=A, 'T'/'C'=A'
 * @param[in] uplo Triangle stored: 'U'=upper, 'L'=lower
 * @param[in] n Order of A, G, Q (n >= 0)
 * @param[in] a N-by-N matrix A, dimension (lda,n)
 * @param[in] lda Leading dimension of A
 * @param[in,out] g N-by-N symmetric matrix G, dimension (ldg,n)
 * @param[in] ldg Leading dimension of G
 * @param[in,out] q N-by-N symmetric matrix Q, dimension (ldq,n)
 * @param[in] ldq Leading dimension of Q
 * @param[out] s 2N-by-2N Hamiltonian/symplectic matrix, dimension (lds,2*n)
 * @param[in] lds Leading dimension of S
 * @param[out] iwork Integer workspace (2*n for discrete)
 * @param[out] dwork Double workspace (6*n for discrete)
 * @param[in] ldwork Workspace size
 * @param[out] info 0=success, 1..n=singular A, n+1=numerically singular A
 */
void sb02ru(
    const char* dico,
    const char* hinv,
    const char* trana,
    const char* uplo,
    const i32 n,
    f64* a,
    const i32 lda,
    f64* g,
    const i32 ldg,
    f64* q,
    const i32 ldq,
    f64* s,
    const i32 lds,
    i32* iwork,
    f64* dwork,
    const i32 ldwork,
    i32* info);

/**
 * @brief Solve linear equations with LU factorization and iterative refinement.
 *
 * Solves op(A)*X = B using LU factorization with optional equilibration
 * and iterative refinement.
 *
 * @param[in] fact Factorization: 'F'=factored, 'N'=factor A, 'E'=equilibrate+factor
 * @param[in] trans Form: 'N'=A*X=B, 'T'/'C'=A'*X=B
 * @param[in] n Order of A (n >= 0)
 * @param[in] nrhs Number of right-hand sides (nrhs >= 0)
 * @param[in,out] a N-by-N matrix A, dimension (lda,n)
 * @param[in] lda Leading dimension of A
 * @param[in,out] af N-by-N LU factors, dimension (ldaf,n)
 * @param[in] ldaf Leading dimension of AF
 * @param[in,out] ipiv Pivot indices, dimension (n)
 * @param[in,out] equed Equilibration type: 'N','R','C','B'
 * @param[in,out] r Row scale factors, dimension (n)
 * @param[in,out] c Column scale factors, dimension (n)
 * @param[in,out] b N-by-NRHS right-hand side B, dimension (ldb,nrhs)
 * @param[in] ldb Leading dimension of B
 * @param[out] x N-by-NRHS solution X, dimension (ldx,nrhs)
 * @param[in] ldx Leading dimension of X
 * @param[out] rcond Reciprocal condition number estimate
 * @param[out] ferr Forward error bounds, dimension (nrhs)
 * @param[out] berr Backward error bounds, dimension (nrhs)
 * @param[out] iwork Integer workspace, dimension (n)
 * @param[out] dwork Double workspace, dimension (4*n)
 * @param[out] info 0=success, <0=invalid param, >0=singular
 */
void mb02pd(
    const char* fact,
    const char* trans,
    const i32 n,
    const i32 nrhs,
    f64* a,
    const i32 lda,
    f64* af,
    const i32 ldaf,
    i32* ipiv,
    char* equed,
    f64* r,
    f64* c,
    f64* b,
    const i32 ldb,
    f64* x,
    const i32 ldx,
    f64* rcond,
    f64* ferr,
    f64* berr,
    i32* iwork,
    f64* dwork,
    i32* info);

/**
 * @brief Solve Lyapunov equation for Cholesky factor of solution.
 *
 * Solves for X = op(U)'*op(U) either the stable continuous-time Lyapunov equation:
 *
 *     op(A)'*X + X*op(A) = -scale^2*op(B)'*op(B)   (DICO='C')
 *
 * or the convergent discrete-time Lyapunov equation:
 *
 *     op(A)'*X*op(A) - X = -scale^2*op(B)'*op(B)   (DICO='D')
 *
 * where op(K) = K or K' (transpose), A is N-by-N, op(B) is M-by-N, and U is
 * the upper triangular Cholesky factor of the solution X. Scale is set <= 1
 * to avoid overflow in X.
 *
 * For continuous-time (DICO='C'): A must be stable (all eigenvalues have
 * negative real parts).
 * For discrete-time (DICO='D'): A must be convergent (all eigenvalues
 * inside the unit circle).
 *
 * Based on Bartels-Stewart method [1] finding Cholesky factor directly without
 * forming the normal matrix op(B)'*op(B) [2].
 *
 * References:
 * [1] Bartels, Stewart. Solution of A'X + XB = C. CACM 15, 820-826, 1972.
 * [2] Hammarling. Numerical solution of stable Lyapunov equation. IMA J. Num. Anal. 2, 303-325, 1982.
 *
 * @param[in] dico Equation type:
 *                 'C' = continuous-time
 *                 'D' = discrete-time
 * @param[in] fact Schur factorization option:
 *                 'F' = A and Q contain Schur factorization (provided by user)
 *                 'N' = Schur factorization will be computed
 * @param[in] trans Form of op(K):
 *                  'N' = op(K) = K (no transpose)
 *                  'T' = op(K) = K' (transpose)
 * @param[in] n Order of matrix A (n >= 0)
 * @param[in] m Number of rows of op(B) (m >= 0)
 * @param[in,out] a N-by-N matrix A, dimension (lda,n)
 *                  If FACT='F': upper quasi-triangular Schur form
 *                  If FACT='N': general matrix, overwritten with Schur form
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in,out] q N-by-N orthogonal matrix Q, dimension (ldq,n)
 *                  If FACT='F': orthogonal matrix from Schur factorization
 *                  If FACT='N': output orthogonal matrix Q
 * @param[in] ldq Leading dimension of Q (ldq >= max(1,n))
 * @param[in,out] b Coefficient matrix B:
 *                  If TRANS='N': M-by-N on entry, dimension (ldb,n)
 *                  If TRANS='T': N-by-M on entry, dimension (ldb,max(m,n))
 *                  On exit: N-by-N upper triangular Cholesky factor U
 * @param[in] ldb Leading dimension of B
 *                If TRANS='N': ldb >= max(1,n,m)
 *                If TRANS='T': ldb >= max(1,n)
 * @param[out] scale Scale factor (0 < scale <= 1) to prevent overflow
 * @param[out] wr Real parts of eigenvalues of A, dimension (n)
 * @param[out] wi Imaginary parts of eigenvalues of A, dimension (n)
 * @param[out] dwork Workspace array, dimension (ldwork)
 *                   On exit: dwork[0] = optimal ldwork
 * @param[in] ldwork Workspace size:
 *                   If m > 0: ldwork >= max(1,4*n)
 *                   If m = 0: ldwork >= 1
 *                   If ldwork = -1: workspace query
 * @param[out] info Exit code:
 *                  0 = success
 *                  1 = nearly singular (warning): perturbed values used
 *                  2 = A not stable/convergent (FACT='N')
 *                  3 = Schur form S not stable/convergent (FACT='F')
 *                  4 = S has >2x2 diagonal block (FACT='F')
 *                  5 = S has 2x2 block with real eigenvalues (FACT='F')
 *                  6 = DGEES failed to converge (FACT='N')
 */
void sb03od(
    const char* dico,
    const char* fact,
    const char* trans,
    const i32 n,
    const i32 m,
    f64* a,
    const i32 lda,
    f64* q,
    const i32 ldq,
    f64* b,
    const i32 ldb,
    f64* scale,
    f64* wr,
    f64* wi,
    f64* dwork,
    const i32 ldwork,
    i32* info);

/**
 * @brief Construct complex plane rotation for Lyapunov solver.
 *
 * Constructs a complex plane rotation such that, for a complex number a and
 * a real number b:
 *
 *     ( conjg(c)   s ) * ( a ) = ( d )
 *     (    -s      c )   ( b )   ( 0 )
 *
 * where d is always real and is overwritten on a, so that on return the
 * imaginary part of a is zero. b is unaltered.
 *
 * @param[in,out] a DOUBLE PRECISION array, dimension (2)
 *                  On entry: a[0] and a[1] are real and imaginary parts of a
 *                  On exit: a[0] contains real part of d, a[1] is set to zero
 * @param[in] b The real number b
 * @param[in] small A small real number. If norm d of [a; b] is smaller than
 *                  small, then the rotation is taken as unit matrix, and
 *                  a[0] and a[1] are set to d and 0, respectively.
 * @param[out] c DOUBLE PRECISION array, dimension (2)
 *               c[0] and c[1] are real and imaginary parts of complex cosine
 * @param[out] s The real sine of the plane rotation
 */
void sb03ov(f64* a, const f64 b, const f64 small, f64* c, f64* s);

/**
 * @brief Solve 2x2 Lyapunov equation for Cholesky factor.
 *
 * Solves for the Cholesky factor U of X, where op(U)'*op(U) = X, either
 * the continuous-time two-by-two Lyapunov equation:
 *
 *     op(S)'*X + X*op(S) = -ISGN*scale^2*op(R)'*op(R)     (DISCR=false)
 *
 * or the discrete-time two-by-two Lyapunov equation:
 *
 *     op(S)'*X*op(S) - X = -ISGN*scale^2*op(R)'*op(R)     (DISCR=true)
 *
 * where op(K) = K or K', S is 2x2 with complex conjugate eigenvalues,
 * R is 2x2 upper triangular, ISGN = -1 or 1, and scale is an output
 * scale factor set <= 1 to avoid overflow in X.
 *
 * Also computes matrices B and A so that:
 *   - B*U = U*S and A*U = scale^2*R (if LTRANS=false), or
 *   - U*B = S*U and U*A = scale^2*R (if LTRANS=true)
 *
 * For continuous-time (DISCR=false), ISGN*S must be stable (eigenvalues
 * have strictly negative real parts). For discrete-time (DISCR=true),
 * if ISGN=1, S must be convergent (eigenvalue moduli < 1); if ISGN=-1,
 * S must be completely divergent (eigenvalue moduli > 1).
 *
 * @param[in] discr Equation type:
 *                  false = continuous-time
 *                  true = discrete-time
 * @param[in] ltrans Form of op(K):
 *                   false = op(K) = K (no transpose)
 *                   true = op(K) = K' (transpose)
 * @param[in] isgn Sign of equation: +1 or -1
 * @param[in,out] s DOUBLE PRECISION array, dimension (lds,2)
 *                  On entry: 2x2 matrix S
 *                  On exit: 2x2 matrix B such that B*U = U*S (LTRANS=false)
 *                           or U*B = S*U (LTRANS=true)
 * @param[in] lds Leading dimension of S (lds >= 2)
 * @param[in,out] r DOUBLE PRECISION array, dimension (ldr,2)
 *                  On entry: 2x2 upper triangular matrix R (R[1,0] not referenced)
 *                  On exit: 2x2 upper triangular Cholesky factor U
 * @param[in] ldr Leading dimension of R (ldr >= 2)
 * @param[out] a DOUBLE PRECISION array, dimension (lda,2)
 *               2x2 upper triangular matrix A satisfying:
 *               A*U/scale = scale*R (LTRANS=false) or
 *               U*A/scale = scale*R (LTRANS=true)
 * @param[in] lda Leading dimension of A (lda >= 2)
 * @param[out] scale Scale factor (0 < scale <= 1) to prevent overflow
 * @param[out] info Exit code:
 *                  0 = success
 *                  1 = near-singular (warning): perturbed values used
 *                  2 = stability requirement not satisfied
 *                  4 = S has real eigenvalues (requires complex conjugate)
 *
 * @note In the interests of speed, this routine does not check all inputs
 *       for errors.
 */
void sb03oy(
    const bool discr,
    const bool ltrans,
    const i32 isgn,
    f64* s, const i32 lds,
    f64* r, const i32 ldr,
    f64* a, const i32 lda,
    f64* scale,
    i32* info);

/**
 * @brief Solve real quasi-triangular Sylvester equation.
 *
 * Solves for the N-by-M matrix X (M = 1 or 2) in:
 *
 *     op(S)'*X + X*op(A) = scale*C   (DISCR = false, continuous)
 *     op(S)'*X*op(A) - X = scale*C   (DISCR = true, discrete)
 *
 * where op(K) = K or K' (transpose), S is an N-by-N block upper triangular
 * matrix with 1x1 and 2x2 blocks on the diagonal (real Schur form), A is an
 * M-by-M matrix. The solution X overwrites C. Scale is set <= 1 to avoid
 * overflow in X.
 *
 * This is a service routine for the Lyapunov solver SB03OT.
 *
 * @param[in] discr Equation type:
 *                  false = continuous: op(S)'*X + X*op(A) = scale*C
 *                  true = discrete: op(S)'*X*op(A) - X = scale*C
 * @param[in] ltrans Form of op(K):
 *                   false = op(K) = K (no transpose)
 *                   true = op(K) = K' (transpose)
 * @param[in] n Order of matrix S, number of rows of X and C (n >= 0)
 * @param[in] m Order of matrix A, number of columns of X and C (m = 1 or 2)
 * @param[in] s N-by-N block upper triangular matrix (real Schur form),
 *              dimension (lds,n). Elements below upper Hessenberg not referenced.
 * @param[in] lds Leading dimension of S (lds >= max(1,n))
 * @param[in] a M-by-M matrix A, dimension (lda,m)
 * @param[in] lda Leading dimension of A (lda >= m)
 * @param[in,out] c N-by-M matrix, dimension (ldc,m)
 *                  On entry: right-hand side matrix C
 *                  On exit: solution matrix X
 * @param[in] ldc Leading dimension of C (ldc >= max(1,n))
 * @param[out] scale Scale factor (0 < scale <= 1) to prevent overflow
 * @param[out] info Exit code:
 *                  0 = success
 *                  1 = S and -A have common eigenvalues (continuous), or
 *                      S and A have eigenvalues with product = 1 (discrete);
 *                      solution computed using slightly perturbed values
 */
void sb03or(
    const bool discr,
    const bool ltrans,
    const i32 n,
    const i32 m,
    const f64* s,
    const i32 lds,
    const f64* a,
    const i32 lda,
    f64* c,
    const i32 ldc,
    f64* scale,
    i32* info);

/**
 * @brief Solve reduced Lyapunov equation for triangular factors.
 *
 * Solves for X = op(U)'*op(U) either the stable continuous-time Lyapunov equation:
 *
 *     op(S)'*X + X*op(S) = -scale^2*op(R)'*op(R)   (continuous)
 *
 * or the convergent discrete-time Lyapunov equation:
 *
 *     op(S)'*X*op(S) - X = -scale^2*op(R)'*op(R)   (discrete)
 *
 * where op(K) = K or K' (transpose), S is an N-by-N block upper triangular matrix
 * with 1x1 or 2x2 blocks on the diagonal (real Schur form), R is an N-by-N upper
 * triangular matrix. The output U is upper triangular and overwrites R. Scale is
 * an output scale factor set <= 1 to avoid overflow.
 *
 * For continuous-time: S must be stable (all eigenvalues have negative real parts).
 * For discrete-time: S must be convergent (all eigenvalues inside unit circle).
 *
 * Based on Bartels-Stewart backward substitution [1] finding Cholesky factor directly
 * without forming the normal matrix op(R)'*op(R) [2].
 *
 * References:
 * [1] Bartels, Stewart. Solution of A'X + XB = C. CACM 15, 820-826, 1972.
 * [2] Hammarling. Numerical solution of stable Lyapunov equation. IMA J. Num. Anal. 2, 303-325, 1982.
 *
 * @param[in] discr Equation type:
 *                  false = continuous: op(S)'*X + X*op(S) = -scale^2*op(R)'*op(R)
 *                  true = discrete: op(S)'*X*op(S) - X = -scale^2*op(R)'*op(R)
 * @param[in] ltrans Form of op(K):
 *                   false = op(K) = K (no transpose)
 *                   true = op(K) = K' (transpose)
 * @param[in] n Order of matrices S and R (n >= 0)
 * @param[in] s N-by-N block upper triangular matrix (real Schur form),
 *              dimension (lds,n). Upper Hessenberg part used; subdiagonal
 *              elements define 2x2 blocks (must correspond to complex
 *              conjugate eigenvalue pairs only).
 * @param[in] lds Leading dimension of S (lds >= max(1,n))
 * @param[in,out] r DOUBLE PRECISION array, dimension (ldr,n)
 *                  On entry: N-by-N upper triangular matrix R
 *                  On exit: N-by-N upper triangular Cholesky factor U
 * @param[in] ldr Leading dimension of R (ldr >= max(1,n))
 * @param[out] scale Scale factor (0 < scale <= 1) to prevent overflow
 * @param[out] dwork Workspace array, dimension (4*n)
 * @param[out] info Exit code:
 *                  0 = success
 *                  1 = near-singular (warning): perturbed values used
 *                  2 = S not stable (continuous) or not convergent (discrete)
 *                  3 = S has >2x2 diagonal block (consecutive non-zero subdiagonals)
 *                  4 = 2x2 block has real eigenvalues (requires complex conjugate)
 */
void sb03ot(
    const bool discr,
    const bool ltrans,
    const i32 n,
    f64* s,
    const i32 lds,
    f64* r,
    const i32 ldr,
    f64* scale,
    f64* dwork,
    i32* info);

/**
 * @brief Solve small (N1xN2, 1<=N1,N2<=2) Sylvester equation.
 *
 * Solves for the N1-by-N2 matrix X (1 <= N1,N2 <= 2) in:
 *
 *     op(TL)*X*op(TR) + ISGN*X = SCALE*B
 *
 * where TL is N1-by-N1, TR is N2-by-N2, B is N1-by-N2, ISGN = 1 or -1,
 * and op(T) = T or T' (transpose).
 *
 * Uses Gaussian elimination with complete pivoting.
 *
 * @param[in] ltranl If true, use TL', otherwise use TL
 * @param[in] ltranr If true, use TR', otherwise use TR
 * @param[in] isgn Sign of equation: +1 or -1
 * @param[in] n1 Order of matrix TL (0, 1, or 2)
 * @param[in] n2 Order of matrix TR (0, 1, or 2)
 * @param[in] tl N1-by-N1 matrix TL, dimension (ldtl,n1)
 * @param[in] ldtl Leading dimension of TL (ldtl >= max(1,n1))
 * @param[in] tr N2-by-N2 matrix TR, dimension (ldtr,n2)
 * @param[in] ldtr Leading dimension of TR (ldtr >= max(1,n2))
 * @param[in] b N1-by-N2 right-hand side matrix B, dimension (ldb,n2)
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n1))
 * @param[out] scale Scale factor (0 < scale <= 1) to prevent overflow
 * @param[out] x N1-by-N2 solution matrix X, dimension (ldx,n2)
 *               Note: X may be identified with B in the calling statement
 * @param[in] ldx Leading dimension of X (ldx >= max(1,n1))
 * @param[out] xnorm Infinity-norm of the solution X
 * @param[out] info Exit code:
 *                  0 = success
 *                  1 = TL and -ISGN*TR have almost reciprocal eigenvalues,
 *                      so TL or TR is perturbed to get nonsingular equation
 *
 * @note This routine does not check inputs for errors (for speed).
 */
void sb04px(
    const bool ltranl,
    const bool ltranr,
    const i32 isgn,
    const i32 n1,
    const i32 n2,
    const f64* tl,
    const i32 ldtl,
    const f64* tr,
    const i32 ldtr,
    const f64* b,
    const i32 ldb,
    f64* scale,
    f64* x,
    const i32 ldx,
    f64* xnorm,
    i32* info);

/**
 * @brief Balance system matrices (A,B,C) using diagonal similarity transformation.
 *
 * Reduces the 1-norm of the system matrix S = [A B; C 0] by balancing.
 * Applies diagonal similarity transformation inv(D)*A*D iteratively
 * to make rows and columns of diag(D,I)^{-1} * S * diag(D,I) as close
 * in norm as possible.
 *
 * The balancing can be performed on:
 *   - S = [A B; C 0] (JOB='A')
 *   - S = [A B]      (JOB='B')
 *   - S = [A; C]     (JOB='C')
 *   - S = A          (JOB='N')
 *
 * @param[in] job Specifies which matrices are involved:
 *                'A' = All matrices (A, B, C)
 *                'B' = B and A matrices only
 *                'C' = C and A matrices only
 *                'N' = A matrix only (B and C not involved)
 * @param[in] n Order of A, rows of B, columns of C (n >= 0)
 * @param[in] m Number of columns of B (m >= 0)
 * @param[in] p Number of rows of C (p >= 0)
 * @param[in,out] maxred On entry: maximum allowed reduction in 1-norm if zero
 *                       rows/columns encountered. If maxred > 0, must be > 1.
 *                       If maxred <= 0, default value 10.0 is used.
 *                       On exit: ratio of original to balanced matrix 1-norm.
 * @param[in,out] a State matrix, dimension (lda,n)
 *                  In: N-by-N matrix A
 *                  Out: Balanced matrix inv(D)*A*D
 * @param[in] lda Leading dimension of A (lda >= max(1,n))
 * @param[in,out] b Input matrix, dimension (ldb,m)
 *                  In: N-by-M matrix B (if m > 0)
 *                  Out: Balanced matrix inv(D)*B
 * @param[in] ldb Leading dimension of B (ldb >= max(1,n) if m > 0, else >= 1)
 * @param[in,out] c Output matrix, dimension (ldc,n)
 *                  In: P-by-N matrix C (if p > 0)
 *                  Out: Balanced matrix C*D
 * @param[in] ldc Leading dimension of C (ldc >= max(1,p))
 * @param[out] scale Scaling factors, dimension (n)
 *                   scale[j] = D(j,j) for j = 0,...,n-1
 * @param[out] info Exit code:
 *                  0 = success
 *                  -i = i-th parameter had illegal value
 */
void tb01id(const char* job, i32 n, i32 m, i32 p, f64* maxred,
            f64* a, i32 lda, f64* b, i32 ldb, f64* c, i32 ldc,
            f64* scale, i32* info);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_H */
