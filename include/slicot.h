#ifndef SLICOT_H
#define SLICOT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Main header for SLICOT C library */

void mb01qd(char type, int32_t m, int32_t n, int32_t kl, int32_t ku,
            double cfrom, double cto, int32_t nbl, const int32_t* nrows,
            double* a, int32_t lda, int32_t* info);

#ifdef __cplusplus
}
#endif

#endif /* SLICOT_H */
