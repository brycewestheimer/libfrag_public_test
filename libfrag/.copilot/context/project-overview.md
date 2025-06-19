# libfrag Project Context

## Project Overview
- **Name**: libfrag
- **Type**: Quantum chemistry library
- **Language**: C++17 with separate header (.hpp) and implementation (.cpp) files
- **Build System**: CMake 3.15+
- **Python Bindings**: pybind11
- **Numerical Library**: xtensor
- **Testing**: doctest (C++), pytest (Python)
- **Documentation**: Sphinx + Doxygen

## Architecture
- **Library structure**: Headers and implementation files separated
- Modern CMake target-based configuration
- Cross-platform support (Linux, macOS, Windows)
- Python package with C++ extension

## Key Dependencies
- xtensor/xtl for numerical computations
- pybind11 for Python bindings
- CMake 3.15+ for build system
- NumPy (Python side)

## Development Standards
- C++17 standard
- const-correctness required
- RAII for resource management
- Exception-safe code
- Comprehensive testing
- Doxygen documentation for C++
- NumPy-style docstrings for Python
- **Clear, descriptive variable names** (no abbreviations)
- **Comprehensive inline comments** explaining algorithms
- **Separate .hpp/.cpp files** for implementation
- **Detailed function documentation** with examples

## Domain-Specific Notes
- Fragment-based quantum chemistry methods
- Molecular system handling
- Numerical stability important
- Performance-critical computations
- Scientific computing best practices
