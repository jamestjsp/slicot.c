#ifndef SLICOT_H
#define SLICOT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Main header for SLICOT C library */

/**
 * @brief Multiply matrix by scalar CTO/CFROM without overflow.
 *
 * Multiplies M-by-N real matrix A by scalar CTO/CFROM, avoiding overflow/underflow
 * as long as the final result can be represented. Supports full, triangular,
 * Hessenberg, and banded matrix types with optional block structure.
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
 * @param[in] cfrom Denominator scalar (must be nonzero)
 * @param[in] cto Numerator scalar
 * @param[in] nbl Number of diagonal blocks (0 = no block structure)
 * @param[in] nrows Block sizes, dimension max(1,nbl)
 * @param[in,out] a Matrix array, dimension (lda,n), column-major
 * @param[in] lda Leading dimension (lda >= max(1,m))
 * @param[out] info Exit code (0 = success, not used in this implementation)
 */
void mb01qd(char type, int32_t m, int32_t n, int32_t kl, int32_t ku,
            double cfrom, double cto, int32_t nbl, const int32_t* nrows,
            double* a, int32_t lda, int32_t* info);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_H */
