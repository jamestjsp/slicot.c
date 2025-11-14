#!/bin/bash
# Setup Python virtual environment for slicot.c development

set -e

VENV_DIR="venv"
PYTHON_CMD="python3"

echo "===================================="
echo "Setting up Python virtual environment"
echo "===================================="

# Check if Python 3 is available
if ! command -v $PYTHON_CMD &> /dev/null; then
    echo "Error: python3 not found. Please install Python 3.8+."
    exit 1
fi

PYTHON_VERSION=$($PYTHON_CMD --version | cut -d' ' -f2)
echo "Using Python version: $PYTHON_VERSION"

# Create virtual environment
if [ -d "$VENV_DIR" ]; then
    echo "Virtual environment already exists at $VENV_DIR"
    read -p "Remove and recreate? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing virtual environment..."
        rm -rf "$VENV_DIR"
    else
        echo "Keeping existing virtual environment."
        exit 0
    fi
fi

echo "Creating virtual environment in $VENV_DIR..."
$PYTHON_CMD -m venv "$VENV_DIR"

# Activate virtual environment
echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install dependencies
echo "Installing dependencies from requirements.txt..."
pip install -r requirements.txt

echo ""
echo "===================================="
echo "Virtual environment setup complete!"
echo "===================================="
echo ""
echo "To activate the virtual environment, run:"
echo "  source venv/bin/activate"
echo ""
echo "To deactivate, run:"
echo "  deactivate"
echo ""
echo "Next steps:"
echo "  1. Activate: source venv/bin/activate"
echo "  2. Build: cmake --preset macos-arm64-debug && cmake --build --preset macos-arm64-debug-build"
echo "  3. Install: pip install -e ."
echo "  4. Test: pytest tests/python/ -v"
echo ""
