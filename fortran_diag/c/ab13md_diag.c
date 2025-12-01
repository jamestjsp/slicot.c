// SPDX-License-Identifier: BSD-3-Clause
//
// AB13MD DIAGNOSTIC PROGRAM - C VERSION
//
// Instrumented version that prints all intermediate values
// for comparison with Fortran implementation.

#include "slicot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#define NMAX 10
#define MMAX 10
#define LDZ NMAX
#define LIWORK (4*MMAX - 2 > NMAX ? 4*MMAX - 2 : NMAX)
#define LDWORK (2*NMAX*NMAX*MMAX - NMAX*NMAX + 9*MMAX*MMAX + NMAX*MMAX + 11*NMAX + 33*MMAX - 11)
#define LZWORK (6*NMAX*NMAX*MMAX + 12*NMAX*NMAX + 6*MMAX + 6*NMAX - 3)

int main(void) {
    int n, m, info;
    int nblock[MMAX], itype[MMAX], iwork[LIWORK];
    double bound;
    double d[NMAX], g[NMAX], x[2*MMAX-1];
    double dwork[LDWORK];
    double complex z[LDZ*NMAX], zwork[LZWORK];

    printf("=== AB13MD C DIAGNOSTIC ===\n\n");

    // Skip header line
    char line[256];
    if (!fgets(line, sizeof(line), stdin)) {
        fprintf(stderr, "ERROR: Failed to read header line\n");
        return 1;
    }

    // Read N, M
    if (scanf("%d %d", &n, &m) != 2) {
        fprintf(stderr, "ERROR: Failed to read N, M\n");
        return 1;
    }

    printf("N = %3d\n", n);
    printf("M = %3d\n", m);

    if (n < 0 || n > NMAX) {
        fprintf(stderr, "ERROR: N out of range: %d\n", n);
        return 1;
    }
    if (m < 0 || m > MMAX) {
        fprintf(stderr, "ERROR: M out of range: %d\n", m);
        return 1;
    }

    // Read NBLOCK
    for (int i = 0; i < m; i++) {
        if (scanf("%d", &nblock[i]) != 1) {
            fprintf(stderr, "ERROR: Failed to read NBLOCK[%d]\n", i);
            return 1;
        }
    }

    // Read ITYPE
    for (int i = 0; i < m; i++) {
        if (scanf("%d", &itype[i]) != 1) {
            fprintf(stderr, "ERROR: Failed to read ITYPE[%d]\n", i);
            return 1;
        }
    }

    // Read Z matrix (complex) - Fortran format: (real,imag) pairs
    // Input is row-wise in the data file
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double re, im;
            // Try to parse Fortran complex format: (re,im) with D notation
            char buf[64];
            if (scanf("%63s", buf) != 1) {
                fprintf(stderr, "ERROR: Failed to read Z(%d,%d)\n", i, j);
                return 1;
            }
            // Replace D with E for C parsing
            for (char* p = buf; *p; p++) {
                if (*p == 'D' || *p == 'd') *p = 'E';
            }
            // Parse (re,im) format
            if (sscanf(buf, "(%lf,%lf)", &re, &im) != 2) {
                fprintf(stderr, "ERROR: Failed to parse complex Z(%d,%d): %s\n", i, j, buf);
                return 1;
            }
            z[i + j*LDZ] = re + im*I;
        }
    }

    printf("\nNBLOCK:\n");
    for (int i = 0; i < m; i++) {
        printf("%5d", nblock[i]);
    }
    printf("\n");

    printf("\nITYPE:\n");
    for (int i = 0; i < m; i++) {
        printf("%5d", itype[i]);
    }
    printf("\n");

    printf("\nINPUT Z matrix (real parts):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf(" %16.9E", creal(z[i + j*LDZ]));
        }
        printf("\n");
    }

    printf("\nINPUT Z matrix (imag parts):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf(" %16.9E", cimag(z[i + j*LDZ]));
        }
        printf("\n");
    }

    // Initialize outputs
    memset(d, 0, sizeof(d));
    memset(g, 0, sizeof(g));
    memset(x, 0, sizeof(x));

    // Call the routine
    int ldz = LDZ;
    int ldwork = LDWORK;
    int lzwork = LZWORK;
    info = ab13md('N', n, z, ldz, m, nblock, itype, x, &bound, d, g,
                  iwork, dwork, ldwork, zwork, lzwork);

    printf("\nINFO = %5d\n", info);

    if (info == 0) {
        printf("\nBOUND:\n");
        printf("%23.16E\n", bound);

        printf("\nD vector:\n");
        for (int i = 0; i < n; i++) {
            printf(" %16.9E", d[i]);
            if ((i + 1) % 5 == 0) printf("\n");
        }
        if (n % 5 != 0) printf("\n");

        printf("\nG vector:\n");
        for (int i = 0; i < n; i++) {
            printf(" %16.9E", g[i]);
            if ((i + 1) % 5 == 0) printf("\n");
        }
        if (n % 5 != 0) printf("\n");

        // Count real blocks (mr)
        int mr = 0;
        for (int i = 0; i < m; i++) {
            if (itype[i] == 1) mr++;
        }
        int mt = m + mr - 1;

        printf("\nX vector (scaling):\n");
        for (int i = 0; i < mt; i++) {
            printf(" %16.9E", x[i]);
            if ((i + 1) % 5 == 0) printf("\n");
        }
        if (mt % 5 != 0) printf("\n");
    } else {
        printf("ERROR: INFO = %3d\n", info);
    }

    return 0;
}
