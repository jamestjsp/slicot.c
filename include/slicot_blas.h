#ifndef SLICOT_BLAS_H
#define SLICOT_BLAS_H

#include "slicot_config.h"
#include "slicot_types.h"

/**
 * @file slicot_blas.h
 * @brief BLAS/LAPACK function declarations with portable symbol naming
 *
 * This header provides a unified interface to BLAS/LAPACK routines that
 * automatically handles different symbol naming conventions:
 * - Lowercase with underscore: dgemm_ (most common)
 * - Lowercase only: dgemm (some systems)
 * - Uppercase: DGEMM (rare)
 *
 * The SLC_FC_* macro is defined by CMake during configuration based on
 * probing the BLAS/LAPACK library.
 *
 * Based on SLICUTLET's approach to BLAS/LAPACK integration.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Symbol mangling macros - defined by CMake configuration */
#ifndef SLC_FC_FUNC
  #if defined(SLC_FC_LOWER_US)
    #define SLC_FC_FUNC(lc, UC) lc##_
  #elif defined(SLC_FC_LOWER)
    #define SLC_FC_FUNC(lc, UC) lc
  #elif defined(SLC_FC_UPPER)
    #define SLC_FC_FUNC(lc, UC) UC
  #else
    /* Default: lowercase with underscore (most common) */
    #define SLC_FC_FUNC(lc, UC) lc##_
  #endif
#endif

/* Integer type for BLAS/LAPACK - future ILP64 support */
#ifdef SLICOT_ILP64
  typedef i64 sl_int;
#else
  typedef i32 sl_int;
#endif

/*
 * BLAS/LAPACK function prototypes
 * Using the "evil hack" from SLICUTLET: temporarily redefine 'int' to 'sl_int'
 * for cleaner prototypes that match Fortran INTEGER parameters.
 */
#define int sl_int

/* BLAS Level 1 - Vector operations */
void SLC_FC_FUNC(dcopy, DCOPY)(const int* n, const f64* x, const int* incx,
                                f64* y, const int* incy);

void SLC_FC_FUNC(daxpy, DAXPY)(const int* n, const f64* alpha, const f64* x,
                                const int* incx, f64* y, const int* incy);

void SLC_FC_FUNC(dscal, DSCAL)(const int* n, const f64* alpha, f64* x,
                                const int* incx);

void SLC_FC_FUNC(dswap, DSWAP)(const int* n, f64* x, const int* incx,
                                f64* y, const int* incy);

f64 SLC_FC_FUNC(ddot, DDOT)(const int* n, const f64* x, const int* incx,
                             const f64* y, const int* incy);

f64 SLC_FC_FUNC(dnrm2, DNRM2)(const int* n, const f64* x, const int* incx);

int SLC_FC_FUNC(idamax, IDAMAX)(const int* n, const f64* x, const int* incx);

/* BLAS Level 2 - Matrix-vector operations */
void SLC_FC_FUNC(dgemv, DGEMV)(const char* trans, const int* m, const int* n,
                                const f64* alpha, const f64* a, const int* lda,
                                const f64* x, const int* incx, const f64* beta,
                                f64* y, const int* incy);

void SLC_FC_FUNC(dtrmv, DTRMV)(const char* uplo, const char* trans, const char* diag,
                                const int* n, const f64* a, const int* lda,
                                f64* x, const int* incx);

/* BLAS Level 3 - Matrix-matrix operations */
void SLC_FC_FUNC(dgemm, DGEMM)(const char* transa, const char* transb,
                                const int* m, const int* n, const int* k,
                                const f64* alpha, const f64* a, const int* lda,
                                const f64* b, const int* ldb, const f64* beta,
                                f64* c, const int* ldc);

void SLC_FC_FUNC(dsyrk, DSYRK)(const char* uplo, const char* trans,
                                const int* n, const int* k, const f64* alpha,
                                const f64* a, const int* lda, const f64* beta,
                                f64* c, const int* ldc);

void SLC_FC_FUNC(dsyr2k, DSYR2K)(const char* uplo, const char* trans,
                                  const int* n, const int* k, const f64* alpha,
                                  const f64* a, const int* lda, const f64* b,
                                  const int* ldb, const f64* beta, f64* c,
                                  const int* ldc);

/* LAPACK - Utilities */
f64 SLC_FC_FUNC(dlamch, DLAMCH)(const char* cmach);

void SLC_FC_FUNC(dlacpy, DLACPY)(const char* uplo, const int* m, const int* n,
                                  const f64* a, const int* lda, f64* b,
                                  const int* ldb);

void SLC_FC_FUNC(dlaset, DLASET)(const char* uplo, const int* m, const int* n,
                                  const f64* alpha, const f64* beta, f64* a,
                                  const int* lda);

/* LAPACK - Factorizations */
void SLC_FC_FUNC(dgetrf, DGETRF)(const int* m, const int* n, f64* a,
                                  const int* lda, int* ipiv, int* info);

void SLC_FC_FUNC(dgeqrf, DGEQRF)(const int* m, const int* n, f64* a,
                                  const int* lda, f64* tau, f64* work,
                                  const int* lwork, int* info);

void SLC_FC_FUNC(dlarfg, DLARFG)(const int* n, f64* alpha, f64* x,
                                  const int* incx, f64* tau);

void SLC_FC_FUNC(dlarf, DLARF)(const char* side, const int* m, const int* n,
                                const f64* v, const int* incv, const f64* tau,
                                f64* c, const int* ldc, f64* work);

/* LAPACK - Singular values and condition estimation */
void SLC_FC_FUNC(dlaic1, DLAIC1)(const int* job, const int* j, const f64* x,
                                  const f64* sest, const f64* w, const f64* gamma,
                                  f64* sestpr, f64* s, f64* c);

f64 SLC_FC_FUNC(dlange, DLANGE)(const char* norm, const int* m, const int* n,
                                  const f64* a, const int* lda, f64* work);

/* LAPACK - Orthogonal transformations */
void SLC_FC_FUNC(dormqr, DORMQR)(const char* side, const char* trans,
                                  const int* m, const int* n, const int* k,
                                  const f64* a, const int* lda, const f64* tau,
                                  f64* c, const int* ldc, f64* work,
                                  const int* lwork, int* info);

void SLC_FC_FUNC(dormrz, DORMRZ)(const char* side, const char* trans,
                                  const int* m, const int* n, const int* k,
                                  const int* l, const f64* a, const int* lda,
                                  const f64* tau, f64* c, const int* ldc,
                                  f64* work, const int* lwork, int* info);

void SLC_FC_FUNC(dtzrzf, DTZRZF)(const int* m, const int* n, f64* a,
                                  const int* lda, f64* tau, f64* work,
                                  const int* lwork, int* info);

#undef int

/* Convenience macros for calling BLAS/LAPACK */
#define SLC_DCOPY    SLC_FC_FUNC(dcopy, DCOPY)
#define SLC_DAXPY    SLC_FC_FUNC(daxpy, DAXPY)
#define SLC_DSCAL    SLC_FC_FUNC(dscal, DSCAL)
#define SLC_DSWAP    SLC_FC_FUNC(dswap, DSWAP)
#define SLC_DDOT     SLC_FC_FUNC(ddot, DDOT)
#define SLC_DNRM2    SLC_FC_FUNC(dnrm2, DNRM2)
#define SLC_IDAMAX   SLC_FC_FUNC(idamax, IDAMAX)
#define SLC_DGEMV    SLC_FC_FUNC(dgemv, DGEMV)
#define SLC_DTRMV    SLC_FC_FUNC(dtrmv, DTRMV)
#define SLC_DGEMM    SLC_FC_FUNC(dgemm, DGEMM)
#define SLC_DSYRK    SLC_FC_FUNC(dsyrk, DSYRK)
#define SLC_DSYR2K   SLC_FC_FUNC(dsyr2k, DSYR2K)
#define SLC_DLAMCH   SLC_FC_FUNC(dlamch, DLAMCH)
#define SLC_DLACPY   SLC_FC_FUNC(dlacpy, DLACPY)
#define SLC_DLASET   SLC_FC_FUNC(dlaset, DLASET)
#define SLC_DGETRF   SLC_FC_FUNC(dgetrf, DGETRF)
#define SLC_DGEQRF   SLC_FC_FUNC(dgeqrf, DGEQRF)
#define SLC_DLARFG   SLC_FC_FUNC(dlarfg, DLARFG)
#define SLC_DLARF    SLC_FC_FUNC(dlarf, DLARF)
#define SLC_DLAIC1   SLC_FC_FUNC(dlaic1, DLAIC1)
#define SLC_DLANGE   SLC_FC_FUNC(dlange, DLANGE)
#define SLC_DORMQR   SLC_FC_FUNC(dormqr, DORMQR)
#define SLC_DORMRZ   SLC_FC_FUNC(dormrz, DORMRZ)
#define SLC_DTZRZF   SLC_FC_FUNC(dtzrzf, DTZRZF)

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_BLAS_H */
