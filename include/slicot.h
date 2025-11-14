#ifndef SLICOT_H
#define SLICOT_H

#include <stdint.h>

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
void mb01qd(char type, int32_t m, int32_t n, int32_t kl, int32_t ku,
            double cfrom, double cto, int32_t nbl, const int32_t* nrows,
            double* a, int32_t lda, int32_t* info);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_H */
