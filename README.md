# slicot.c

C11 translation of SLICOT (Subroutine Library In Control Theory).

## Build

Using CMake presets:

```bash
# Configure
cmake --preset macos-arm64-debug

# Build
cmake --build --preset macos-arm64-debug-build

# Test
ctest --preset macos-arm64-debug-test
```

Available presets:
- `macos-x64-debug` / `macos-x64-release`
- `macos-arm64-debug` / `macos-arm64-release`

## Options

- `BUILD_SHARED_LIBS`: Build shared libraries (default: ON)
- `BUILD_TESTING`: Build tests (default: ON)
- `BUILD_EXAMPLES`: Build examples (default: ON)
