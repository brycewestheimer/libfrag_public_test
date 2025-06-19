# Agent-Based Development Guidelines for libfrag

This document outlines best practices for developing libfrag using AI coding assistants like GitHub Copilot, Claude, ChatGPT, and other AI tools.

## Project Context

**libfrag** is a quantum chemistry library for fragment-based calculations with:
- Modern C++17 core implementation with separate header and source files
- Python bindings via pybind11
- CMake build system
- Cross-platform support
- Scientific computing focus with xtensor integration
- Emphasis on clear, self-documenting code with meaningful variable names

## Development Workflow with AI Assistants

### 1. Code Generation Guidelines

#### For C++ Development:
- Always specify C++17 features when requesting code
- Emphasize const-correctness and RAII principles
- Request modern CMake patterns (target-based, not directory-based)
- Ask for exception-safe code with proper error handling
- Prefer standard library algorithms over raw loops
- Request Doxygen-style documentation for public APIs
- **Use clear, descriptive variable names** that explain purpose/meaning
- **Add comprehensive comments** explaining algorithms and business logic
- **Separate headers and implementation files** (.hpp/.cpp structure)
- **Document all function parameters** and return values
- **Explain complex algorithms** with inline comments
- **Use meaningful class and function names** that describe their purpose

#### For Python Bindings:
- Use pybind11 modern syntax (py::module_ instead of py::module)
- Request NumPy integration examples when dealing with arrays
- Ask for proper docstrings following NumPy style
- Ensure type hints are included in Python code
- Request error handling that maps C++ exceptions to Python exceptions
- **Use descriptive parameter names** in Python function signatures
- **Include comprehensive docstrings** with examples and parameter descriptions
- **Add type hints** for all function parameters and return values
- **Document scientific/mathematical concepts** in docstrings

#### For CMake:
- Use modern CMake 3.15+ features
- Request target-based configurations
- Ask for proper find_package usage
- Ensure cross-platform compatibility

### 2. Code Review Prompts

When reviewing AI-generated code, ask the assistant to check for:

```
Please review this code for:
1. Memory safety and RAII compliance
2. const-correctness
3. Exception safety
4. Performance considerations
5. Quantum chemistry domain appropriateness
6. Cross-platform compatibility
7. CMake best practices (if applicable)
8. **Variable naming clarity and meaningfulness**
9. **Code documentation and comment quality**
10. **Separation of interface and implementation**
11. **Function and class naming descriptiveness**
```

### 3. Testing Strategy

When requesting test code:
- Ask for both unit tests (doctest/gtest) and integration tests
- Request property-based testing for numerical methods
- Ask for Python pytest fixtures for complex setups
- Request benchmark tests for performance-critical code
- Ask for cross-platform test considerations

### 4. Documentation Requests

- Request both inline documentation and user guides
- Ask for mathematical notation in docstrings where appropriate
- Request examples that demonstrate quantum chemistry concepts
- Ask for API documentation that explains both C++ and Python interfaces

### 5. Domain-Specific Considerations

When working on quantum chemistry features:
- Provide context about molecular systems and fragment-based methods
- Explain any specific quantum chemistry terminology
- Ask for numerically stable implementations
- Request validation against known chemical systems
- Ask for proper handling of molecular symmetry and periodicity

## File Organization

- `chats/` - Conversation logs with AI assistants
- `prompts/` - Reusable prompt templates
- `guidelines/` - This and other development guidelines
- `context/` - Project context files for sharing with assistants
- `templates/` - Code templates for common patterns

## Prompt Templates

### New Feature Development
```
I'm working on libfrag, a C++17 quantum chemistry library for fragment-based calculations. 
I need to implement [feature description]. The code should:
- Follow modern C++17 practices
- Be compatible with our CMake 3.15+ build system
- Include proper error handling
- Have Python bindings via pybind11
- Include unit tests
- Be documented with Doxygen comments

[Specific requirements...]
```

### Code Review Request
```
Please review this C++17 code for libfrag (quantum chemistry library):
- Check for memory safety and const-correctness
- Verify exception safety
- Assess performance for numerical computations
- Ensure compatibility with our pybind11 Python bindings
- Validate CMake integration

[Code to review...]
```

## Code Quality and Documentation Standards

### Variable Naming Conventions
- Use **descriptive, full words** instead of abbreviations
- Examples:
  - ✅ `molecular_fragment_count` instead of ❌ `mfc` or `n`
  - ✅ `orbital_coefficient_matrix` instead of ❌ `coeff` or `c`
  - ✅ `total_energy_hartrees` instead of ❌ `energy` or `e`
  - ✅ `convergence_threshold` instead of ❌ `tol` or `eps`

### Function and Class Naming
- Use **verb-noun combinations** for functions: `calculate_fragment_energy()`, `optimize_geometry()`
- Use **descriptive class names**: `MolecularFragment`, `QuantumChemicalCalculator`
- Avoid generic names like `Manager`, `Handler`, `Processor`

### Documentation Requirements

#### C++ Documentation (Doxygen):
```cpp
/**
 * @brief Calculates the total electronic energy of a molecular fragment
 * 
 * This function performs a self-consistent field calculation to determine
 * the ground state electronic energy of the given molecular fragment using
 * the specified basis set and exchange-correlation functional.
 * 
 * @param molecular_fragment The fragment containing atomic coordinates and charges
 * @param basis_set_name Name of the basis set (e.g., "cc-pVDZ", "6-31G*")
 * @param xcfunctional_type Exchange-correlation functional type
 * @param convergence_threshold SCF convergence criterion in hartrees
 * 
 * @return Total electronic energy in atomic units (hartrees)
 * 
 * @throws std::invalid_argument if basis set is not recognized
 * @throws std::runtime_error if SCF fails to converge
 * 
 * @note This function assumes the fragment geometry is already optimized
 * @warning Large basis sets may require significant memory (>4GB for >100 atoms)
 * 
 * @example
 * ```cpp
 * MolecularFragment water_fragment = create_water_molecule();
 * double energy = calculate_fragment_energy(water_fragment, "cc-pVDZ", 
 *                                         XCFunctional::B3LYP, 1e-8);
 * ```
 */
double calculate_fragment_energy(const MolecularFragment& molecular_fragment,
                               const std::string& basis_set_name,
                               XCFunctional xcfunctional_type,
                               double convergence_threshold = 1e-6);
```

#### Python Documentation (NumPy Style):
```python
def calculate_fragment_energy(molecular_fragment: MolecularFragment,
                            basis_set_name: str,
                            xcfunctional_type: str = "B3LYP",
                            convergence_threshold: float = 1e-6) -> float:
    """Calculate the total electronic energy of a molecular fragment.
    
    Performs a self-consistent field calculation to determine the ground state
    electronic energy using the specified quantum chemical method.
    
    Parameters
    ----------
    molecular_fragment : MolecularFragment
        Fragment containing atomic coordinates, charges, and multiplicity
    basis_set_name : str
        Name of the basis set (e.g., 'cc-pVDZ', '6-31G*')
    xcfunctional_type : str, default 'B3LYP'
        Exchange-correlation functional name
    convergence_threshold : float, default 1e-6
        SCF convergence criterion in hartrees
        
    Returns
    -------
    float
        Total electronic energy in atomic units (hartrees)
        
    Raises
    ------
    ValueError
        If basis set name is not recognized
    RuntimeError
        If SCF calculation fails to converge
        
    Notes
    -----
    This function assumes the molecular geometry is already optimized.
    Large basis sets may require significant computational resources.
    
    The SCF algorithm uses DIIS acceleration for faster convergence.
    
    Examples
    --------
    >>> water = create_water_molecule()
    >>> energy = calculate_fragment_energy(water, 'cc-pVDZ')
    >>> print(f"Water energy: {energy:.6f} hartrees")
    Water energy: -76.241234 hartrees
    
    >>> # Using different functional and tighter convergence
    >>> energy_pbe = calculate_fragment_energy(water, '6-31G*', 'PBE', 1e-8)
    """
```

### Implementation File Organization
- **Header files (.hpp)**: Class declarations, inline functions, templates
- **Source files (.cpp)**: Function implementations, static data
- **Separate compilation units** for major classes/functionality
- **Clear separation** between public interface and private implementation

## Best Practices

1. **Always provide project context** when starting new conversations
2. **Be specific about C++17 and CMake versions** in requests
3. **Mention quantum chemistry domain** when relevant
4. **Request cross-platform considerations** for all code
5. **Ask for both C++ and Python examples** when implementing new features
6. **Request proper error handling** for numerical edge cases
7. **Ask for performance considerations** in tight loops
8. **Request comprehensive documentation** that explains both implementation and usage
9. **Emphasize clear variable naming** - avoid abbreviations and single letters
10. **Ask for detailed comments** explaining algorithms and business logic
11. **Request separate .hpp/.cpp files** for classes and complex functions
12. **Ask for inline documentation** of mathematical formulas and scientific concepts
13. **Request examples in docstrings** showing typical usage patterns

## Version Control

- Keep conversation logs for major features
- Tag important prompts and responses
- Document successful prompt patterns
- Note any AI-generated code that needed significant modifications
