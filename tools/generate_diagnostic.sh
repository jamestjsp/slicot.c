#!/usr/bin/env bash
# SPDX-License-Identifier: BSD-3-Clause
#
# Generate diagnostic program skeleton for a SLICOT routine
#
# Usage: ./tools/generate_diagnostic.sh ROUTINE_NAME
#
# Example: ./tools/generate_diagnostic.sh AB01MD
#
# This script:
# 1. Creates fortran_diag/fortran/<routine>_diag.f from template
# 2. Creates fortran_diag/c/<routine>_diag.c from template
# 3. Updates fortran_diag/CMakeLists.txt to add the routine
# 4. Prints next steps for manual completion

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check arguments
if [ $# -ne 1 ]; then
    echo -e "${RED}ERROR: Missing routine name${NC}"
    echo "Usage: $0 ROUTINE_NAME"
    echo "Example: $0 AB01MD"
    exit 1
fi

ROUTINE=$(echo "$1" | tr '[:lower:]' '[:upper:]')
routine_lower=$(echo "$1" | tr '[:upper:]' '[:lower:]')

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Paths
FORTRAN_DIAG_DIR="$PROJECT_ROOT/fortran_diag"
FORTRAN_TEMPLATE="$FORTRAN_DIAG_DIR/templates/skeleton.f"
C_TEMPLATE="$FORTRAN_DIAG_DIR/templates/skeleton.c"
FORTRAN_OUT="$FORTRAN_DIAG_DIR/fortran/${routine_lower}_diag.f"
C_OUT="$FORTRAN_DIAG_DIR/c/${routine_lower}_diag.c"
CMAKE_FILE="$FORTRAN_DIAG_DIR/CMakeLists.txt"

# SLICOT-Reference paths
SLICOT_REF="$PROJECT_ROOT/SLICOT-Reference"
EXAMPLE_PROGRAM="$SLICOT_REF/examples/T${ROUTINE}.f"
DATA_FILE="$SLICOT_REF/examples/data/${ROUTINE}.dat"
EXPECTED_FILE="$SLICOT_REF/examples/results/${ROUTINE}.res"

echo -e "${GREEN}=== Generating diagnostic for $ROUTINE ===${NC}\n"

# Check if templates exist
if [ ! -f "$FORTRAN_TEMPLATE" ]; then
    echo -e "${RED}ERROR: Fortran template not found: $FORTRAN_TEMPLATE${NC}"
    exit 1
fi

if [ ! -f "$C_TEMPLATE" ]; then
    echo -e "${RED}ERROR: C template not found: $C_TEMPLATE${NC}"
    exit 1
fi

# Check if SLICOT-Reference files exist
echo "Checking SLICOT-Reference files..."
REF_FILES_FOUND=true

if [ ! -f "$EXAMPLE_PROGRAM" ]; then
    echo -e "${YELLOW}WARNING: Reference example not found: $EXAMPLE_PROGRAM${NC}"
    REF_FILES_FOUND=false
fi

if [ ! -f "$DATA_FILE" ]; then
    echo -e "${YELLOW}WARNING: Test data file not found: $DATA_FILE${NC}"
    REF_FILES_FOUND=false
fi

if [ ! -f "$EXPECTED_FILE" ]; then
    echo -e "${YELLOW}WARNING: Expected results not found: $EXPECTED_FILE${NC}"
    REF_FILES_FOUND=false
fi

# Check if output files already exist
if [ -f "$FORTRAN_OUT" ]; then
    echo -e "${YELLOW}WARNING: $FORTRAN_OUT already exists${NC}"
    read -p "Overwrite? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping Fortran file generation"
        FORTRAN_OUT=""
    fi
fi

if [ -f "$C_OUT" ]; then
    echo -e "${YELLOW}WARNING: $C_OUT already exists${NC}"
    read -p "Overwrite? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping C file generation"
        C_OUT=""
    fi
fi

# Generate Fortran diagnostic
if [ -n "$FORTRAN_OUT" ]; then
    echo "Generating $FORTRAN_OUT..."
    sed "s/{{ROUTINE}}/$ROUTINE/g" "$FORTRAN_TEMPLATE" > "$FORTRAN_OUT"
    echo -e "${GREEN}✓ Created $FORTRAN_OUT${NC}"
fi

# Generate C diagnostic
if [ -n "$C_OUT" ]; then
    echo "Generating $C_OUT..."
    sed -e "s/{{ROUTINE}}/$ROUTINE/g" \
        -e "s/{{routine_lower}}/$routine_lower/g" \
        "$C_TEMPLATE" > "$C_OUT"
    echo -e "${GREEN}✓ Created $C_OUT${NC}"
fi

# Update CMakeLists.txt
echo ""
echo "Updating $CMAKE_FILE..."

# Check if routine is already in CMakeLists.txt
if grep -q "add_slicot_diagnostic($ROUTINE)" "$CMAKE_FILE"; then
    echo -e "${YELLOW}✓ $ROUTINE already in CMakeLists.txt (no changes needed)${NC}"
else
    # Find the line with "# Add diagnostics for implemented routines" and add after it
    # Use sed to insert after the last add_slicot_diagnostic line
    if grep -q "add_slicot_diagnostic" "$CMAKE_FILE"; then
        # Insert after the last add_slicot_diagnostic line
        sed -i.bak "/add_slicot_diagnostic/a\\
add_slicot_diagnostic($ROUTINE)
" "$CMAKE_FILE"
        rm "$CMAKE_FILE.bak"
        echo -e "${GREEN}✓ Added $ROUTINE to CMakeLists.txt${NC}"
    else
        echo -e "${YELLOW}WARNING: Could not automatically update CMakeLists.txt${NC}"
        echo "Please manually add: add_slicot_diagnostic($ROUTINE)"
    fi
fi

# Print next steps
echo ""
echo -e "${GREEN}=== NEXT STEPS ===${NC}"
echo ""
echo "1. ${YELLOW}Implement parsing logic:${NC}"

if [ "$REF_FILES_FOUND" = true ]; then
    echo "   Reference example: $EXAMPLE_PROGRAM"
    echo "   Look for READ statements in the example program"
    echo ""
    echo "   Edit $FORTRAN_OUT"
    echo "   Edit $C_OUT"
else
    echo -e "   ${RED}WARNING: Reference files not found - you'll need to determine the data format manually${NC}"
fi

echo ""
echo "2. ${YELLOW}Build and test:${NC}"
echo "   cmake --preset macos-arm64-debug -DBUILD_FORTRAN_DIAG=ON"
echo "   cmake --build --preset macos-arm64-debug-build --target ${routine_lower}_diag_all"
echo ""
echo "3. ${YELLOW}View results:${NC}"
echo "   cat build/macos-arm64-debug/fortran_diag/${routine_lower}_fortran.txt"
echo "   cat build/macos-arm64-debug/fortran_diag/${routine_lower}_c.txt"
echo "   cat build/macos-arm64-debug/fortran_diag/output/${routine_lower}_diff.txt"
echo ""
echo -e "${GREEN}Done!${NC}"
