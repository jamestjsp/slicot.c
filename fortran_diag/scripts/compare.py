#!/usr/bin/env python3
"""
Compare SLICOT Fortran and C diagnostic outputs.

Parses trace files from instrumented programs and compares
numerical values to identify where implementations diverge.
"""

import sys
import re
import numpy as np
from typing import Dict, List, Tuple


def parse_output(filename: str) -> Dict[str, np.ndarray]:
    """
    Parse diagnostic output file and extract matrices/vectors.

    Args:
        filename: Path to trace file

    Returns:
        Dictionary mapping matrix/vector names to numpy arrays
    """
    data = {}
    current_name = None
    current_data = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Detect matrix/vector headers (lines ending with ':')
            if line.endswith(':'):
                # Save previous matrix/vector if exists
                if current_name and current_data:
                    data[current_name] = np.array(current_data)

                # Start new matrix/vector
                current_name = line.rstrip(':')
                current_data = []
                continue

            # Try to parse numerical data
            # Match scientific notation: [-+]?d.dddE[-+]dd
            numbers = re.findall(r'[-+]?\d+\.\d+[Ee][-+]?\d+', line)
            if numbers:
                row = [float(x) for x in numbers]
                current_data.append(row)

        # Save last matrix/vector
        if current_name and current_data:
            data[current_name] = np.array(current_data)

    return data


def extract_scalars(filename: str) -> Dict[str, float]:
    """
    Extract scalar values (INFO, SCALE, N, M) from output.

    Args:
        filename: Path to trace file

    Returns:
        Dictionary mapping scalar names to values
    """
    scalars = {}

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # Match patterns like "INFO = 0" or "SCALE = 1.23E+00"
            match = re.match(r'(INFO|SCALE|N|M)\s*=\s*(.+)', line)
            if match:
                name = match.group(1)
                value_str = match.group(2).strip()

                # Try to parse as number
                try:
                    if 'E' in value_str or 'e' in value_str:
                        value = float(value_str)
                    else:
                        value = int(value_str)
                    scalars[name] = value
                except ValueError:
                    pass

    return scalars


def compare_values(fortran_data: Dict[str, np.ndarray],
                   c_data: Dict[str, np.ndarray],
                   fortran_scalars: Dict[str, float],
                   c_scalars: Dict[str, float],
                   tolerance: float = 1e-12) -> Tuple[bool, List[str]]:
    """
    Compare Fortran and C outputs and report differences.

    Args:
        fortran_data: Fortran matrices/vectors
        c_data: C matrices/vectors
        fortran_scalars: Fortran scalar values
        c_scalars: C scalar values
        tolerance: Absolute tolerance for numerical comparison

    Returns:
        Tuple of (all_match, differences)
        - all_match: True if everything matches within tolerance
        - differences: List of difference descriptions
    """
    differences = []
    all_match = True

    print("=" * 80)
    print("FORTRAN vs C DIAGNOSTIC COMPARISON")
    print("=" * 80)

    # Compare scalars
    print("\n" + "=" * 80)
    print("SCALAR VALUES")
    print("=" * 80)

    all_scalar_keys = set(fortran_scalars.keys()) | set(c_scalars.keys())

    for key in sorted(all_scalar_keys):
        if key not in fortran_scalars:
            msg = f"  {key}: MISSING in Fortran"
            print(msg)
            differences.append(msg)
            all_match = False
        elif key not in c_scalars:
            msg = f"  {key}: MISSING in C"
            print(msg)
            differences.append(msg)
            all_match = False
        else:
            f_val = fortran_scalars[key]
            c_val = c_scalars[key]

            if isinstance(f_val, int) and isinstance(c_val, int):
                # Integer comparison
                if f_val == c_val:
                    print(f"  {key}: MATCH ({f_val})")
                else:
                    msg = f"  {key}: MISMATCH (Fortran={f_val}, C={c_val})"
                    print(msg)
                    differences.append(msg)
                    all_match = False
            else:
                # Floating point comparison
                diff = abs(f_val - c_val)
                if diff <= tolerance:
                    print(f"  {key}: MATCH ({f_val:.6e})")
                else:
                    msg = f"  {key}: MISMATCH (Fortran={f_val:.6e}, C={c_val:.6e}, diff={diff:.6e})"
                    print(msg)
                    differences.append(msg)
                    all_match = False

    # Compare matrices/vectors
    print("\n" + "=" * 80)
    print("MATRIX/VECTOR VALUES")
    print("=" * 80)

    all_keys = set(fortran_data.keys()) | set(c_data.keys())

    for key in sorted(all_keys):
        print(f"\n{key}:")
        print("-" * 80)

        if key not in fortran_data:
            msg = f"  MISSING in Fortran output"
            print(msg)
            differences.append(f"{key}: {msg}")
            all_match = False
            continue

        if key not in c_data:
            msg = f"  MISSING in C output"
            print(msg)
            differences.append(f"{key}: {msg}")
            all_match = False
            continue

        f_array = fortran_data[key]
        c_array = c_data[key]

        # Check shapes
        if f_array.shape != c_array.shape:
            msg = f"  SHAPE MISMATCH: Fortran {f_array.shape} vs C {c_array.shape}"
            print(msg)
            differences.append(f"{key}: {msg}")
            all_match = False
            continue

        # Flatten for vector output
        if len(f_array.shape) == 1:
            f_flat = f_array
            c_flat = c_array
        else:
            f_flat = f_array.flatten()
            c_flat = c_array.flatten()

        # Compute differences
        diff = np.abs(f_flat - c_flat)
        max_diff = np.max(diff)
        rms_diff = np.sqrt(np.mean(diff**2))
        max_idx = np.argmax(diff)

        if max_diff > tolerance:
            msg = f"  MISMATCH DETECTED!"
            print(msg)
            print(f"  Max absolute difference: {max_diff:.6e}")
            print(f"  RMS difference: {rms_diff:.6e}")
            print(f"  Location of max diff: index {max_idx}")
            print(f"    Fortran value: {f_flat[max_idx]:.16e}")
            print(f"    C value: {c_flat[max_idx]:.16e}")

            # Show first few mismatches
            mismatches = np.where(diff > tolerance)[0]
            if len(mismatches) > 0:
                print(f"  Total mismatches: {len(mismatches)}/{len(diff)}")
                print(f"  First 5 mismatch indices: {mismatches[:5].tolist()}")

            differences.append(f"{key}: max_diff={max_diff:.6e}, rms={rms_diff:.6e}")
            all_match = False
        else:
            print(f"  MATCH (max diff: {max_diff:.6e})")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    if all_match:
        print("✓ All values MATCH within tolerance")
    else:
        print(f"✗ Found {len(differences)} MISMATCHES")
        print("\nDifferences:")
        for diff in differences:
            print(f"  - {diff}")

    print("=" * 80)

    return all_match, differences


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: compare.py <fortran_output> <c_output>", file=sys.stderr)
        print("\nCompares diagnostic traces from Fortran and C implementations.", file=sys.stderr)
        print("Reports numerical differences to identify divergence points.", file=sys.stderr)
        sys.exit(1)

    fortran_file = sys.argv[1]
    c_file = sys.argv[2]

    print(f"Comparing:")
    print(f"  Fortran: {fortran_file}")
    print(f"  C:       {c_file}")
    print()

    try:
        # Parse both files
        fortran_data = parse_output(fortran_file)
        c_data = parse_output(c_file)
        fortran_scalars = extract_scalars(fortran_file)
        c_scalars = extract_scalars(c_file)

        # Compare
        all_match, differences = compare_values(
            fortran_data, c_data,
            fortran_scalars, c_scalars
        )

        # Exit code: 0 if match, 1 if mismatch
        sys.exit(0 if all_match else 1)

    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(2)


if __name__ == '__main__':
    main()
