/*
 * SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 1996-2025, The SLICOT Team (original Fortran77 code)
 * Copyright (c) 2025, slicot.c contributors (C11 translation)
 */

#ifndef SLICOT_BLAS_H
#define SLICOT_BLAS_H

#include "slicot_config.h"
#include "slicot_types.h"
#include <stdbool.h>

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

f64 SLC_FC_FUNC(dasum, DASUM)(const int* n, const f64* x, const int* incx);

/* BLAS Level 2 - Matrix-vector operations */
void SLC_FC_FUNC(dgemv, DGEMV)(const char* trans, const int* m, const int* n,
                                const f64* alpha, const f64* a, const int* lda,
                                const f64* x, const int* incx, const f64* beta,
                                f64* y, const int* incy);

void SLC_FC_FUNC(dtrmv, DTRMV)(const char* uplo, const char* trans, const char* diag,
                                const int* n, const f64* a, const int* lda,
                                f64* x, const int* incx);

void SLC_FC_FUNC(dtrsv, DTRSV)(const char* uplo, const char* trans, const char* diag,
                                const int* n, const f64* a, const int* lda,
                                f64* x, const int* incx);

void SLC_FC_FUNC(dger, DGER)(const int* m, const int* n, const f64* alpha,
                              const f64* x, const int* incx, const f64* y,
                              const int* incy, f64* a, const int* lda);

/* BLAS Level 3 - Matrix-matrix operations */
void SLC_FC_FUNC(dtrmm, DTRMM)(const char* side, const char* uplo, const char* transa,
                                const char* diag, const int* m, const int* n,
                                const f64* alpha, const f64* a, const int* lda,
                                f64* b, const int* ldb);

void SLC_FC_FUNC(dtrsm, DTRSM)(const char* side, const char* uplo, const char* transa,
                                const char* diag, const int* m, const int* n,
                                const f64* alpha, const f64* a, const int* lda,
                                f64* b, const int* ldb);

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

void SLC_FC_FUNC(dsymm, DSYMM)(const char* side, const char* uplo,
                                const int* m, const int* n, const f64* alpha,
                                const f64* a, const int* lda, const f64* b,
                                const int* ldb, const f64* beta, f64* c,
                                const int* ldc);

/* LAPACK - Utilities */
f64 SLC_FC_FUNC(dlamch, DLAMCH)(const char* cmach);

void SLC_FC_FUNC(dlacpy, DLACPY)(const char* uplo, const int* m, const int* n,
                                  const f64* a, const int* lda, f64* b,
                                  const int* ldb);

void SLC_FC_FUNC(dlapmt, DLAPMT)(const int* forwrd, const int* m, const int* n,
                                  f64* x, const int* ldx, int* k);

void SLC_FC_FUNC(dlarnv, DLARNV)(const int* idist, int* seed, const int* n,
                                  f64* x);

void SLC_FC_FUNC(dlaset, DLASET)(const char* uplo, const int* m, const int* n,
                                  const f64* alpha, const f64* beta, f64* a,
                                  const int* lda);

/* LAPACK - Factorizations */
void SLC_FC_FUNC(dgetrf, DGETRF)(const int* m, const int* n, f64* a,
                                  const int* lda, int* ipiv, int* info);

void SLC_FC_FUNC(dgetrs, DGETRS)(const char* trans, const int* n, const int* nrhs,
                                  const f64* a, const int* lda, const int* ipiv,
                                  f64* b, const int* ldb, int* info);

void SLC_FC_FUNC(dgetri, DGETRI)(const int* n, f64* a, const int* lda,
                                  const int* ipiv, f64* work, const int* lwork,
                                  int* info);

void SLC_FC_FUNC(dgesv, DGESV)(const int* n, const int* nrhs, f64* a,
                                const int* lda, int* ipiv, f64* b,
                                const int* ldb, int* info);

void SLC_FC_FUNC(dgecon, DGECON)(const char* norm, const int* n, const f64* a,
                                  const int* lda, const f64* anorm, f64* rcond,
                                  f64* work, int* iwork, int* info);

void SLC_FC_FUNC(dpotrf, DPOTRF)(const char* uplo, const int* n, f64* a,
                                  const int* lda, int* info);

void SLC_FC_FUNC(dsytrf, DSYTRF)(const char* uplo, const int* n, f64* a,
                                  const int* lda, int* ipiv, f64* work,
                                  const int* lwork, int* info);

void SLC_FC_FUNC(dsysv, DSYSV)(const char* uplo, const int* n, const int* nrhs,
                                f64* a, const int* lda, int* ipiv, f64* b,
                                const int* ldb, f64* work, const int* lwork,
                                int* info);

void SLC_FC_FUNC(dpocon, DPOCON)(const char* uplo, const int* n, const f64* a,
                                  const int* lda, const f64* anorm, f64* rcond,
                                  f64* work, int* iwork, int* info);

void SLC_FC_FUNC(dsycon, DSYCON)(const char* uplo, const int* n, const f64* a,
                                  const int* lda, const int* ipiv, const f64* anorm,
                                  f64* rcond, f64* work, int* iwork, int* info);

void SLC_FC_FUNC(dsytrs, DSYTRS)(const char* uplo, const int* n, const int* nrhs,
                                  const f64* a, const int* lda, const int* ipiv,
                                  f64* b, const int* ldb, int* info);

void SLC_FC_FUNC(dpotrs, DPOTRS)(const char* uplo, const int* n, const int* nrhs,
                                  const f64* a, const int* lda,
                                  f64* b, const int* ldb, int* info);

void SLC_FC_FUNC(dtrcon, DTRCON)(const char* norm, const char* uplo,
                                  const char* diag, const int* n, const f64* a,
                                  const int* lda, f64* rcond, f64* work,
                                  int* iwork, int* info);

f64 SLC_FC_FUNC(dlansy, DLANSY)(const char* norm, const char* uplo, const int* n,
                                 const f64* a, const int* lda, f64* work);

void SLC_FC_FUNC(dgeqrf, DGEQRF)(const int* m, const int* n, f64* a,
                                  const int* lda, f64* tau, f64* work,
                                  const int* lwork, int* info);

void SLC_FC_FUNC(dgeqp3, DGEQP3)(const int* m, const int* n, f64* a,
                                  const int* lda, int* jpvt, f64* tau,
                                  f64* work, const int* lwork, int* info);

void SLC_FC_FUNC(dgerqf, DGERQF)(const int* m, const int* n, f64* a,
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

f64 SLC_FC_FUNC(dlanhs, DLANHS)(const char* norm, const int* n,
                                 const f64* a, const int* lda, f64* work);

f64 SLC_FC_FUNC(dlapy2, DLAPY2)(const f64* x, const f64* y);

f64 SLC_FC_FUNC(dlapy3, DLAPY3)(const f64* x, const f64* y, const f64* z);
f64 SLC_FC_FUNC(dlantr, DLANTR)(const char* norm, const char* uplo, const char* diag,
                                const int* m, const int* n, const f64* a,
                                const int* lda, f64* work);

void SLC_FC_FUNC(dlascl, DLASCL)(const char* type, const int* kl, const int* ku,
                                  const f64* cfrom, const f64* cto, const int* m,
                                  const int* n, f64* a, const int* lda, int* info);

void SLC_FC_FUNC(dlabad, DLABAD)(f64* small, f64* large);

void SLC_FC_FUNC(dlag2, DLAG2)(const f64* a, const int* lda, const f64* b, const int* ldb,
                                const f64* safmin, f64* scale1, f64* scale2,
                                f64* wr1, f64* wr2, f64* wi);

void SLC_FC_FUNC(dlasv2, DLASV2)(const f64* f, const f64* g, const f64* h,
                                  f64* ssmin, f64* ssmax, f64* snr, f64* csr,
                                  f64* snl, f64* csl);

void SLC_FC_FUNC(dladiv, DLADIV)(const f64* a, const f64* b, const f64* c, const f64* d,
                                  f64* p, f64* q);

void SLC_FC_FUNC(dlanv2, DLANV2)(f64* a, f64* b, f64* c, f64* d,
                                  f64* rt1r, f64* rt1i, f64* rt2r, f64* rt2i,
                                  f64* cs, f64* sn);

void SLC_FC_FUNC(zlarfg, ZLARFG)(const int* n, double complex* alpha, double complex* x,
                                  const int* incx, double complex* tau);

void SLC_FC_FUNC(zstein, ZSTEIN)(const int* n, const f64* d, const f64* e, const int* m,
                                  const f64* w, const int* iblock, const int* isplit,
                                  double complex* z, const int* ldz, f64* work,
                                  int* iwork, int* ifail, int* info);

/* Complex BLAS */
void SLC_FC_FUNC(zcopy, ZCOPY)(const int* n, const double complex* x, const int* incx,
                                double complex* y, const int* incy);

void SLC_FC_FUNC(zgemv, ZGEMV)(const char* trans, const int* m, const int* n,
                                const double complex* alpha, const double complex* a,
                                const int* lda, const double complex* x, const int* incx,
                                const double complex* beta, double complex* y, const int* incy);

void SLC_FC_FUNC(zgemm, ZGEMM)(const char* transa, const char* transb,
                                const int* m, const int* n, const int* k,
                                const double complex* alpha, const double complex* a,
                                const int* lda, const double complex* b, const int* ldb,
                                const double complex* beta, double complex* c, const int* ldc);

/* Complex LAPACK */
void SLC_FC_FUNC(zlacpy, ZLACPY)(const char* uplo, const int* m, const int* n,
                                  const double complex* a, const int* lda,
                                  double complex* b, const int* ldb);

void SLC_FC_FUNC(zlascl, ZLASCL)(const char* type, const int* kl, const int* ku,
                                  const f64* cfrom, const f64* cto, const int* m,
                                  const int* n, double complex* a, const int* lda, int* info);

f64 SLC_FC_FUNC(zlange, ZLANGE)(const char* norm, const int* m, const int* n,
                                 const double complex* a, const int* lda, f64* work);

void SLC_FC_FUNC(zgetrf, ZGETRF)(const int* m, const int* n, double complex* a,
                                  const int* lda, int* ipiv, int* info);

void SLC_FC_FUNC(zgetri, ZGETRI)(const int* n, double complex* a, const int* lda,
                                  const int* ipiv, double complex* work, const int* lwork,
                                  int* info);

void SLC_FC_FUNC(zgees, ZGEES)(const char* jobvs, const char* sort,
                                int (*select)(const double complex*),
                                const int* n, double complex* a, const int* lda,
                                int* sdim, double complex* w, double complex* vs,
                                const int* ldvs, double complex* work, const int* lwork,
                                f64* rwork, int* bwork, int* info);

void SLC_FC_FUNC(zgesvd, ZGESVD)(const char* jobu, const char* jobvt,
                                  const int* m, const int* n, double complex* a,
                                  const int* lda, f64* s, double complex* u, const int* ldu,
                                  double complex* vt, const int* ldvt, double complex* work,
                                  const int* lwork, f64* rwork, int* info);

void SLC_FC_FUNC(xerbla, XERBLA)(const char* srname, const int* info);

void SLC_FC_FUNC(dgetc2, DGETC2)(const int* n, f64* a, const int* lda,
                                  int* ipiv, int* jpiv, int* info);

void SLC_FC_FUNC(dgesc2, DGESC2)(const int* n, const f64* a, const int* lda,
                                  f64* rhs, const int* ipiv, const int* jpiv,
                                  f64* scale);

/* LAPACK - Sylvester equations */
void SLC_FC_FUNC(dlasy2, DLASY2)(const int* ltranl, const int* ltranr,
                                  const int* isgn, const int* n1, const int* n2,
                                  const f64* tl, const int* ldtl,
                                  const f64* tr, const int* ldtr,
                                  const f64* b, const int* ldb,
                                  f64* scale, f64* x, const int* ldx,
                                  f64* xnorm, int* info);

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

void SLC_FC_FUNC(dlartg, DLARTG)(const f64* f, const f64* g, f64* cs, f64* sn,
                                  f64* r);

void SLC_FC_FUNC(drot, DROT)(const int* n, f64* dx, const int* incx,
                              f64* dy, const int* incy, const f64* c,
                              const f64* s);

void SLC_FC_FUNC(dsyevx, DSYEVX)(const char* jobz, const char* range, const char* uplo,
                                  const int* n, f64* a, const int* lda,
                                  const f64* vl, const f64* vu, const int* il,
                                  const int* iu, const f64* abstol, int* m,
                                  f64* w, f64* z, const int* ldz, f64* work,
                                  const int* lwork, int* iwork, int* ifail,
                                  int* info);

void SLC_FC_FUNC(dsyev, DSYEV)(const char* jobz, const char* uplo, const int* n,
                                f64* a, const int* lda, f64* w, f64* work,
                                const int* lwork, int* info);

void SLC_FC_FUNC(dgges, DGGES)(const char* jobvsl, const char* jobvsr,
                                const char* sort, int (*selctg)(),
                                const int* n, f64* a, const int* lda,
                                f64* b, const int* ldb, int* sdim,
                                f64* alphar, f64* alphai, f64* beta,
                                f64* vsl, const int* ldvsl, f64* vsr,
                                const int* ldvsr, f64* work, const int* lwork,
                                int* bwork, int* info);

void SLC_FC_FUNC(dgees, DGEES)(const char* jobvs, const char* sort,
                                int (*select)(const f64*, const f64*),
                                const int* n, f64* a, const int* lda,
                                int* sdim, f64* wr, f64* wi,
                                f64* vs, const int* ldvs, f64* work,
                                const int* lwork, int* bwork, int* info);

/* LAPACK - Bidiagonal decomposition and SVD */
void SLC_FC_FUNC(dgebrd, DGEBRD)(const int* m, const int* n, f64* a,
                                  const int* lda, f64* d, f64* e, f64* tauq,
                                  f64* taup, f64* work, const int* lwork,
                                  int* info);

void SLC_FC_FUNC(dorgbr, DORGBR)(const char* vect, const int* m, const int* n,
                                  const int* k, f64* a, const int* lda,
                                  const f64* tau, f64* work, const int* lwork,
                                  int* info);

void SLC_FC_FUNC(dbdsqr, DBDSQR)(const char* uplo, const int* n, const int* ncvt,
                                  const int* nru, const int* ncc, f64* d, f64* e,
                                  f64* vt, const int* ldvt, f64* u, const int* ldu,
                                  f64* c, const int* ldc, f64* work, int* info);

void SLC_FC_FUNC(dgelss, DGELSS)(const int* m, const int* n, const int* nrhs,
                                  f64* a, const int* lda, f64* b, const int* ldb,
                                  f64* s, const f64* rcond, int* rank, f64* work,
                                  const int* lwork, int* info);

void SLC_FC_FUNC(dgelsy, DGELSY)(const int* m, const int* n, const int* nrhs,
                                  f64* a, const int* lda, f64* b, const int* ldb,
                                  int* jpvt, const f64* rcond, int* rank, f64* work,
                                  const int* lwork, int* info);

void SLC_FC_FUNC(dtrtrs, DTRTRS)(const char* uplo, const char* trans, const char* diag,
                                  const int* n, const int* nrhs, const f64* a,
                                  const int* lda, f64* b, const int* ldb, int* info);

void SLC_FC_FUNC(dormrz, DORMRZ)(const char* side, const char* trans,
                                  const int* m, const int* n, const int* k, const int* l,
                                  const f64* a, const int* lda, const f64* tau,
                                  f64* c, const int* ldc, f64* work, const int* lwork,
                                  int* info);

int SLC_FC_FUNC(ilaenv, ILAENV)(const int* ispec, const char* name, const char* opts,
                                const int* n1, const int* n2, const int* n3, const int* n4);

#undef int

/* Convenience macros for calling BLAS/LAPACK */
#define SLC_DCOPY    SLC_FC_FUNC(dcopy, DCOPY)
#define SLC_DAXPY    SLC_FC_FUNC(daxpy, DAXPY)
#define SLC_DSCAL    SLC_FC_FUNC(dscal, DSCAL)
#define SLC_DSWAP    SLC_FC_FUNC(dswap, DSWAP)
#define SLC_DDOT     SLC_FC_FUNC(ddot, DDOT)
#define SLC_DNRM2    SLC_FC_FUNC(dnrm2, DNRM2)
#define SLC_IDAMAX   SLC_FC_FUNC(idamax, IDAMAX)
#define SLC_DASUM    SLC_FC_FUNC(dasum, DASUM)
#define SLC_DGEMV    SLC_FC_FUNC(dgemv, DGEMV)
#define SLC_DTRMV    SLC_FC_FUNC(dtrmv, DTRMV)
#define SLC_DTRSV    SLC_FC_FUNC(dtrsv, DTRSV)
#define SLC_DGER     SLC_FC_FUNC(dger, DGER)
#define SLC_DTRMM    SLC_FC_FUNC(dtrmm, DTRMM)
#define SLC_DTRSM    SLC_FC_FUNC(dtrsm, DTRSM)
#define SLC_DGEMM    SLC_FC_FUNC(dgemm, DGEMM)
#define SLC_DSYRK    SLC_FC_FUNC(dsyrk, DSYRK)
#define SLC_DSYR2K   SLC_FC_FUNC(dsyr2k, DSYR2K)
#define SLC_DSYMM    SLC_FC_FUNC(dsymm, DSYMM)
#define SLC_DLAMCH   SLC_FC_FUNC(dlamch, DLAMCH)
#define SLC_DLACPY   SLC_FC_FUNC(dlacpy, DLACPY)
#define SLC_DLAPMT   SLC_FC_FUNC(dlapmt, DLAPMT)
#define SLC_DLARNV   SLC_FC_FUNC(dlarnv, DLARNV)
#define SLC_DLASET   SLC_FC_FUNC(dlaset, DLASET)
#define SLC_DGETRF   SLC_FC_FUNC(dgetrf, DGETRF)
#define SLC_DGETRS   SLC_FC_FUNC(dgetrs, DGETRS)
#define SLC_DGETRI   SLC_FC_FUNC(dgetri, DGETRI)
#define SLC_DGESV    SLC_FC_FUNC(dgesv, DGESV)
#define SLC_DGECON   SLC_FC_FUNC(dgecon, DGECON)
#define SLC_DPOTRF   SLC_FC_FUNC(dpotrf, DPOTRF)
#define SLC_DSYTRF   SLC_FC_FUNC(dsytrf, DSYTRF)
#define SLC_DSYSV    SLC_FC_FUNC(dsysv, DSYSV)
#define SLC_DPOCON   SLC_FC_FUNC(dpocon, DPOCON)
#define SLC_DSYCON   SLC_FC_FUNC(dsycon, DSYCON)
#define SLC_DSYTRS   SLC_FC_FUNC(dsytrs, DSYTRS)
#define SLC_DLANSY   SLC_FC_FUNC(dlansy, DLANSY)
#define SLC_DGEQRF   SLC_FC_FUNC(dgeqrf, DGEQRF)
#define SLC_DGEQP3   SLC_FC_FUNC(dgeqp3, DGEQP3)
#define SLC_DGERQF   SLC_FC_FUNC(dgerqf, DGERQF)
#define SLC_DLARFG   SLC_FC_FUNC(dlarfg, DLARFG)
#define SLC_DLARF    SLC_FC_FUNC(dlarf, DLARF)
#define SLC_DLAIC1   SLC_FC_FUNC(dlaic1, DLAIC1)
#define SLC_DLANGE   SLC_FC_FUNC(dlange, DLANGE)
#define SLC_DLANHS   SLC_FC_FUNC(dlanhs, DLANHS)
#define SLC_DORMQR   SLC_FC_FUNC(dormqr, DORMQR)
#define SLC_DORMRZ   SLC_FC_FUNC(dormrz, DORMRZ)
#define SLC_DTZRZF   SLC_FC_FUNC(dtzrzf, DTZRZF)
#define SLC_DLAPY2   SLC_FC_FUNC(dlapy2, DLAPY2)
#define SLC_DLAPY3   SLC_FC_FUNC(dlapy3, DLAPY3)
#define SLC_DLASCL   SLC_FC_FUNC(dlascl, DLASCL)
#define SLC_DLABAD   SLC_FC_FUNC(dlabad, DLABAD)
#define SLC_DLAG2    SLC_FC_FUNC(dlag2, DLAG2)
#define SLC_DLASV2   SLC_FC_FUNC(dlasv2, DLASV2)
#define SLC_DLADIV   SLC_FC_FUNC(dladiv, DLADIV)
#define SLC_DLANV2   SLC_FC_FUNC(dlanv2, DLANV2)
#define SLC_ZLARFG   SLC_FC_FUNC(zlarfg, ZLARFG)
#define SLC_ZSTEIN   SLC_FC_FUNC(zstein, ZSTEIN)
#define SLC_ZCOPY    SLC_FC_FUNC(zcopy, ZCOPY)
#define SLC_ZGEMV    SLC_FC_FUNC(zgemv, ZGEMV)
#define SLC_ZGEMM    SLC_FC_FUNC(zgemm, ZGEMM)
#define SLC_ZLACPY   SLC_FC_FUNC(zlacpy, ZLACPY)
#define SLC_ZLASCL   SLC_FC_FUNC(zlascl, ZLASCL)
#define SLC_ZLANGE   SLC_FC_FUNC(zlange, ZLANGE)
#define SLC_ZGETRF   SLC_FC_FUNC(zgetrf, ZGETRF)
#define SLC_ZGETRI   SLC_FC_FUNC(zgetri, ZGETRI)
#define SLC_ZGEES    SLC_FC_FUNC(zgees, ZGEES)
#define SLC_ZGESVD   SLC_FC_FUNC(zgesvd, ZGESVD)
#define SLC_XERBLA   SLC_FC_FUNC(xerbla, XERBLA)
#define SLC_DGETC2   SLC_FC_FUNC(dgetc2, DGETC2)
#define SLC_DGESC2   SLC_FC_FUNC(dgesc2, DGESC2)
#define SLC_DLASY2   SLC_FC_FUNC(dlasy2, DLASY2)
#define SLC_DLARTG   SLC_FC_FUNC(dlartg, DLARTG)
#define SLC_DROT     SLC_FC_FUNC(drot, DROT)
#define SLC_DSYEVX   SLC_FC_FUNC(dsyevx, DSYEVX)
#define SLC_DSYEV    SLC_FC_FUNC(dsyev, DSYEV)
#define SLC_DGGES    SLC_FC_FUNC(dgges, DGGES)
#define SLC_DGEES    SLC_FC_FUNC(dgees, DGEES)
#define SLC_DGEBRD   SLC_FC_FUNC(dgebrd, DGEBRD)
#define SLC_DORGBR   SLC_FC_FUNC(dorgbr, DORGBR)
#define SLC_DBDSQR   SLC_FC_FUNC(dbdsqr, DBDSQR)
#define SLC_DGELSS   SLC_FC_FUNC(dgelss, DGELSS)
#define SLC_DLANTR   SLC_FC_FUNC(dlantr, DLANTR)
#define SLC_DTRCON   SLC_FC_FUNC(dtrcon, DTRCON)
#define SLC_DPOTRS   SLC_FC_FUNC(dpotrs, DPOTRS)
#define SLC_DTRTRS   SLC_FC_FUNC(dtrtrs, DTRTRS)
#define SLC_DGELSY   SLC_FC_FUNC(dgelsy, DGELSY)
#define SLC_DTZRZF   SLC_FC_FUNC(dtzrzf, DTZRZF)
#define SLC_DORMRZ   SLC_FC_FUNC(dormrz, DORMRZ)
#define SLC_ILAENV   SLC_FC_FUNC(ilaenv, ILAENV)

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_BLAS_H */
