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
    void (*fcn)(i32*, i32, i32, const i32*, i32, const f64*, i32, const f64*, i32,
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
    const i32* ipar,
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

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_H */
