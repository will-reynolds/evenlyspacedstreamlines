# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python package for generating evenly-spaced streamlines from an orientation field on a triangulated 3D surface. The core algorithm is implemented in C++ for performance, with Cython bindings and a Python wrapper interface.

## Architecture

- **Core C++ Engine**: Performance-critical algorithms in `evenlyspacedstreamlines/*.{cpp,h}` files
- **Cython Interface**: `runengine.pyx` provides the bridge between C++ and Python
- **Python Wrapper**: `wrapper.py` contains the main API (`evenly_spaced_streamlines` function)
- **Build System**: Uses setuptools with Cython extensions, configured in `setup.py`

The main entry point is the `evenly_spaced_streamlines()` function which takes vertices, triangles, orientation vectors, and a radius parameter to generate streamlines.

## Development Commands

### Building the Package
```bash
# Build wheel package
make all

# Local development build (builds extension in-place)
make local

# Build source distribution
make sdist
```

### Testing
```bash
# Quick functionality test
python -c "from evenlyspacedstreamlines import test; test()"

# Run C++ unit tests
cd tests/unittests && make auto

# Run specific test categories
cd tests/unittests && make triangularmesh
cd tests/unittests && make wrapper

# Run benchmarks
cd tests/benchmarks && make streamlines

# Check for memory leaks
cd tests/benchmarks && make memoryleaks
```

### Development Workflow
1. For local development, use `make local` to build extensions in-place
2. Set `PYTHONPATH=.:$PYTHONPATH` when running tests from subdirectories
3. The package requires numpy and cython as build dependencies
4. Uses OpenMP for parallelization (compile flags: `-fopenmp` on Linux, `/openmp` on Windows)

### Key Files to Understand
- `evenlyspacedstreamlines/wrapper.py`: Main Python API
- `evenlyspacedstreamlines/engine.cpp`: Core streamline generation algorithm  
- `evenlyspacedstreamlines/runengine.pyx`: Cython interface to C++ engine
- `setup.py`: Build configuration with platform-specific compiler flags

The test suite includes unit tests (C++), benchmarks, visual tests, and Python wrapper tests in separate subdirectories under `tests/`.