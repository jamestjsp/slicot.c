#!/usr/bin/env python3
import numpy as np
from slicot import mb04id

# Simple 4x3 matrix with p=1
n, m, p = 4, 3, 1
a = np.array([
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
    [0.0, 10., 11.]
], dtype=float, order='F')

a_orig = a.copy()
print("Original A:")
print(a_orig)
print()

a_out, tau, info = mb04id(n, m, p, a)
print("After MB04ID:")
print("info =", info)
print("tau =", tau)
print("A (contains R and Householder vectors):")
print(a_out)
print()

# For column 0 (i=0 < p=1):
# DLARFG(N-P=3, A(0,0), A(1,0), 1, TAU(0))
# Creates reflector of length 3, stored at rows 0:3, column 0
# The reflector acts on rows 0:3

# For column 1 (i=1 >= p=1):
# DGEQRF handles remaining (n-p) x (m-p) = 3x2 submatrix at A(p:n, p:m)
# Standard QR on A[1:4, 1:3]

# Extract the reflector vectors
print("Reconstructing Q:")
print("=" * 50)

# Column 0 reflector (structured)
print(f"\nColumn 0 (i < p): reflector length = {n-p}")
print(f"  Stored at: A[{0}:{0+(n-p)}, {0}]")
print(f"  Values: diagonal={a_out[0,0]}, below diagonal={a_out[1:0+(n-p), 0]}")
print(f"  tau[0] = {tau[0]}")

# Columns 1-2 reflectors (from DGEQRF on remaining submatrix)
for i in range(1, min(n, m)):
    print(f"\nColumn {i} (i >= p): standard reflector")
    print(f"  Stored at: A[{i}:{n}, {i}]")
    print(f"  Values: diagonal={a_out[i,i]}, below diagonal={a_out[i+1:n, i]}")
    print(f"  tau[{i}] = {tau[i]}")
