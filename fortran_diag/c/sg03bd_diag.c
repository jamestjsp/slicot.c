// SPDX-License-Identifier: BSD-3-Clause
//
// SG03BD DIAGNOSTIC PROGRAM - C VERSION
//
// Instrumented version that prints all intermediate values
// for comparison with Fortran implementation.

#include "slicot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define NMAX 20
#define LDWORK (6*NMAX)

static void print_matrix(const char* name, const double* mat, int rows, int cols, int ld) {
    printf("\n%s:\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %23.16E", mat[i + j*ld]);
        }
        printf("\n");
    }
}

static void print_vector(const char* name, const double* vec, int n) {
    printf("\n%s:\n", name);
    for (int i = 0; i < n; i++) {
        printf(" %23.16E", vec[i]);
    }
    printf("\n");
}

// Helper function: Compute Frobenius norm for validation
static double frobenius_norm(const double* mat, int rows, int cols, int ld) {
    double sum = 0.0;
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            double val = mat[i + j*ld];
            sum += val * val;
        }
    }
    return sqrt(sum);
}

// Helper function: Compute checksum for validation
static double checksum(const double* mat, int rows, int cols, int ld) {
    double sum = 0.0;
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            sum += mat[i + j*ld];
        }
    }
    return sum;
}

// Helper function: Find max absolute element
static double max_abs(const double* mat, int rows, int cols, int ld) {
    double max_val = 0.0;
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            double abs_val = fabs(mat[i + j*ld]);
            if (abs_val > max_val) {
                max_val = abs_val;
            }
        }
    }
    return max_val;
}

int main(void) {
    char dico, fact, trans;
    int n, m, info;
    double scale;

    // Allocate arrays (column-major storage)
    double a[NMAX*NMAX], e[NMAX*NMAX], q[NMAX*NMAX], z[NMAX*NMAX];
    double b[NMAX*NMAX];
    double alphar[NMAX], alphai[NMAX], beta[NMAX];
    double dwork[LDWORK];

    // Initialize arrays
    memset(a, 0, sizeof(a));
    memset(e, 0, sizeof(e));
    memset(q, 0, sizeof(q));
    memset(z, 0, sizeof(z));
    memset(b, 0, sizeof(b));
    memset(alphar, 0, sizeof(alphar));
    memset(alphai, 0, sizeof(alphai));
    memset(beta, 0, sizeof(beta));

    printf("=== SG03BD C DIAGNOSTIC ===\n\n");

    // Skip header line
    char line[256];
    if (!fgets(line, sizeof(line), stdin)) {
        fprintf(stderr, "ERROR: Failed to read header line\n");
        return 1;
    }

    // Read parameters
    if (scanf("%d %d %c %c %c", &n, &m, &dico, &fact, &trans) != 5) {
        fprintf(stderr, "ERROR: Failed to read parameters\n");
        return 1;
    }

    // Convert to uppercase
    dico = toupper((unsigned char)dico);
    fact = toupper((unsigned char)fact);
    trans = toupper((unsigned char)trans);

    printf("N = %d\n", n);
    printf("M = %d\n", m);
    printf("DICO = %c\n", dico);
    printf("FACT = %c\n", fact);
    printf("TRANS = %c\n", trans);

    if (n < 0 || n > NMAX) {
        fprintf(stderr, "ERROR: N out of range: %d\n", n);
        return 1;
    }
    if (m < 0 || m > NMAX) {
        fprintf(stderr, "ERROR: M out of range: %d\n", m);
        return 1;
    }

    // Read A matrix (Fortran row-wise input → C column-major storage)
    // Input is row i, columns j=1..n
    // Store as column-major: a[i + j*NMAX]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (scanf("%lf", &a[i + j*NMAX]) != 1) {
                fprintf(stderr, "ERROR: Failed to read A(%d,%d)\n", i, j);
                return 1;
            }
        }
    }

    // Read E matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (scanf("%lf", &e[i + j*NMAX]) != 1) {
                fprintf(stderr, "ERROR: Failed to read E(%d,%d)\n", i, j);
                return 1;
            }
        }
    }

    // Read Q, Z if FACT='F'
    if (fact == 'F') {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (scanf("%lf", &q[i + j*NMAX]) != 1) {
                    fprintf(stderr, "ERROR: Failed to read Q(%d,%d)\n", i, j);
                    return 1;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (scanf("%lf", &z[i + j*NMAX]) != 1) {
                    fprintf(stderr, "ERROR: Failed to read Z(%d,%d)\n", i, j);
                    return 1;
                }
            }
        }
    }

    // Read B matrix (size depends on TRANS)
    if (trans == 'T') {
        // B is N x M
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (scanf("%lf", &b[i + j*NMAX]) != 1) {
                    fprintf(stderr, "ERROR: Failed to read B(%d,%d)\n", i, j);
                    return 1;
                }
            }
        }
    } else {
        // B is M x N
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (scanf("%lf", &b[i + j*NMAX]) != 1) {
                    fprintf(stderr, "ERROR: Failed to read B(%d,%d)\n", i, j);
                    return 1;
                }
            }
        }
    }

    // Print input matrices
    print_matrix("INPUT A matrix", a, n, n, NMAX);
    print_matrix("INPUT E matrix", e, n, n, NMAX);

    if (fact == 'F') {
        print_matrix("INPUT Q matrix", q, n, n, NMAX);
        print_matrix("INPUT Z matrix", z, n, n, NMAX);
    }

    if (trans == 'T') {
        print_matrix("INPUT B matrix", b, n, m, NMAX);
    } else {
        print_matrix("INPUT B matrix", b, m, n, NMAX);
    }

    // Print validation metrics for inputs
    printf("\n=== INPUT VALIDATION ===\n");
    printf("A Frobenius norm: %23.16E\n", frobenius_norm(a, n, n, NMAX));
    printf("A checksum:       %23.16E\n", checksum(a, n, n, NMAX));
    printf("A max abs:        %23.16E\n", max_abs(a, n, n, NMAX));
    printf("E Frobenius norm: %23.16E\n", frobenius_norm(e, n, n, NMAX));
    printf("E checksum:       %23.16E\n", checksum(e, n, n, NMAX));
    printf("E max abs:        %23.16E\n", max_abs(e, n, n, NMAX));
    if (trans == 'T') {
        printf("B Frobenius norm: %23.16E\n", frobenius_norm(b, n, m, NMAX));
        printf("B checksum:       %23.16E\n", checksum(b, n, m, NMAX));
        printf("B max abs:        %23.16E\n", max_abs(b, n, m, NMAX));
    } else {
        printf("B Frobenius norm: %23.16E\n", frobenius_norm(b, m, n, NMAX));
        printf("B checksum:       %23.16E\n", checksum(b, m, n, NMAX));
        printf("B max abs:        %23.16E\n", max_abs(b, m, n, NMAX));
    }

    // Call sg03bd
    sg03bd(&dico, &fact, &trans, n, m, a, NMAX, e, NMAX, q, NMAX, z, NMAX,
           b, NMAX, &scale, alphar, alphai, beta, dwork, LDWORK, &info);

    // Print intermediate B matrix for debugging
    printf("\nDEBUG: B after sg03bv/before transformations:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf(" %23.16E", b[i + j*NMAX]);
        }
        printf("\n");
    }

    // Print results
    printf("\nINFO = %d\n", info);

    if (info == 0) {
        printf("\nSCALE = %23.16E\n", scale);

        print_matrix("OUTPUT U matrix (Cholesky factor)", b, n, n, NMAX);
        print_vector("ALPHAR eigenvalues", alphar, n);
        print_vector("ALPHAI eigenvalues", alphai, n);
        print_vector("BETA values", beta, n);
        print_matrix("OUTPUT A matrix (Schur form)", a, n, n, NMAX);
        print_matrix("OUTPUT E matrix (triangular)", e, n, n, NMAX);
        print_matrix("OUTPUT Q matrix (orthogonal)", q, n, n, NMAX);
        print_matrix("OUTPUT Z matrix (orthogonal)", z, n, n, NMAX);

        // Print validation metrics for outputs
        printf("\n=== OUTPUT VALIDATION ===\n");
        printf("U Frobenius norm: %23.16E\n", frobenius_norm(b, n, n, NMAX));
        printf("U checksum:       %23.16E\n", checksum(b, n, n, NMAX));
        printf("U max abs:        %23.16E\n", max_abs(b, n, n, NMAX));
        printf("A Frobenius norm: %23.16E\n", frobenius_norm(a, n, n, NMAX));
        printf("A checksum:       %23.16E\n", checksum(a, n, n, NMAX));
        printf("E Frobenius norm: %23.16E\n", frobenius_norm(e, n, n, NMAX));
        printf("E checksum:       %23.16E\n", checksum(e, n, n, NMAX));
        printf("Q Frobenius norm: %23.16E\n", frobenius_norm(q, n, n, NMAX));
        printf("Z Frobenius norm: %23.16E\n", frobenius_norm(z, n, n, NMAX));
    } else {
        printf("SG03BD failed - no output matrices\n");
    }

    return info;
}
