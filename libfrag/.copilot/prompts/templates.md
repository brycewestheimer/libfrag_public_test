# Common Prompt Templates for libfrag Development

## Feature Implementation Template

```
I'm working on libfrag, a C++17 quantum chemistry library for fragment-based calculations.

**Project Context:**
- C++17 with separate .hpp/.cpp files and CMake 3.15+ build system
- Python bindings via pybind11
- Uses xtensor for numerical computations
- Cross-platform support required
- Scientific computing focus
- **Emphasis on clear variable naming and comprehensive documentation**

**Request:**
I need to implement [FEATURE_DESCRIPTION].

**Requirements:**
- Modern C++17 practices with separate .hpp/.cpp files
- Exception-safe and const-correct
- Include unit tests (doctest for C++, pytest for Python)
- Comprehensive Doxygen documentation with examples
- Python bindings if applicable
- CMake integration
- Performance considerations for numerical code
- **Use descriptive, meaningful variable names (no abbreviations)**
- **Add detailed inline comments explaining algorithms and logic**
- **Include comprehensive function documentation with parameter descriptions**

**Specific Details:**
[ADD_SPECIFIC_REQUIREMENTS]
```

## Code Review Template

```
Please review this code for libfrag (C++17 quantum chemistry library):

**Review Criteria:**
1. Memory safety and RAII compliance
2. Const-correctness
3. Exception safety
4. Performance for numerical computations
5. Cross-platform compatibility
6. CMake best practices (if applicable)
7. Python binding efficiency (if applicable)
8. Quantum chemistry domain appropriateness
9. **Variable naming clarity and descriptiveness**
10. **Code documentation quality and completeness**
11. **Proper separation of header and implementation files**
12. **Function and class naming meaningfulness**

**Code:**
[PASTE_CODE_HERE]
```

## CMake Help Template

```
I need CMake help for libfrag (C++17 quantum chemistry library).

**Build System Requirements:**
- CMake 3.15+
- Modern target-based configuration
- Cross-platform support (Linux, macOS, Windows)
- Integration with xtensor, pybind11
- Support for both C++ library and Python extension builds
- Package configuration for find_package() support
- **Separate compilation for .hpp/.cpp file pairs**

**Specific Request:**
[DESCRIBE_CMAKE_NEED]
```

## Python Binding Template

```
I need pybind11 bindings for libfrag (C++17 quantum chemistry library).

**Context:**
- C++17 backend with numerical computations
- xtensor for array operations
- Need NumPy integration
- Scientific computing use case
- Cross-platform Python package

**Requirements:**
- Modern pybind11 syntax
- Proper exception handling (C++ -> Python)
- NumPy array integration
- Type hints where applicable
- NumPy-style docstrings with comprehensive examples
- Performance-efficient bindings
- **Descriptive parameter names in all function signatures**
- **Detailed docstrings explaining scientific concepts**
- **Complete type annotations for all functions**

**Specific Request:**
[DESCRIBE_BINDING_NEED]
```

## Testing Template

```
I need tests for libfrag (C++17 quantum chemistry library).

**Testing Requirements:**
- C++: doctest framework
- Python: pytest framework
- Both unit and integration tests
- Numerical accuracy validation
- Cross-platform test compatibility
- Property-based testing for numerical methods
- Performance benchmarks where appropriate

**Code to Test:**
[DESCRIBE_FUNCTIONALITY]
```

## Documentation Template

```
I need comprehensive documentation for libfrag (C++17 quantum chemistry library).

**Documentation Requirements:**
- Doxygen-style comments for C++ code
- NumPy-style docstrings for Python code
- Clear variable names that explain their purpose
- Inline comments explaining complex algorithms
- Mathematical formulas and scientific concepts explained
- Examples showing typical usage patterns
- Parameter descriptions with units and constraints
- Return value descriptions with expected ranges
- Exception documentation with when they occur
- Performance notes and memory requirements

**Code to Document:**
[PASTE_CODE_HERE]

**Specific Focus:**
[DESCRIBE_WHAT_NEEDS_SPECIAL_ATTENTION]
```
