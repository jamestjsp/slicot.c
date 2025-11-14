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

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_H */
