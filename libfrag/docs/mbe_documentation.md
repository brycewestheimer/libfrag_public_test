# Many-Body Expansion (MBE) Module

The Many-Body Expansion module provides comprehensive functionality for fragment-based quantum chemistry calculations using the MBE method.

## Overview

Many-Body Expansion is a method that expresses the total energy of a molecular system as a sum of contributions from individual fragments and their interactions:

```
E_total = Σ E_i + Σ E_ij + Σ E_ijk + ...
```

Where:
- `E_i` are 1-body (monomer) energies
- `E_ij` are 2-body interaction energies
- `E_ijk` are 3-body interaction energies
- And so on...

## Key Components

### MBEConfig
Configuration class for MBE calculations with support for:
- Fragmentation schemes (atomic, molecular, functional group, custom, distance-based)
- Truncation methods (order-based, distance-based, energy-based)
- QM method specifications
- Convergence criteria

### MBEResults
Results storage and analysis with features for:
- Energy decomposition by N-body order
- Convergence analysis and error estimation
- Performance statistics and timing
- Export to JSON, CSV, and formatted reports

### MBEFragmentGenerator
Fragment combination generation supporting:
- Multiple fragmentation algorithms
- Distance and energy-based filtering
- Combinatorial generation with duplicate removal
- Functional group recognition (placeholder)

### MBECalculator
Main calculation orchestrator providing:
- Parallel quantum chemistry calculation management
- Progress tracking and logging
- Caching for computational efficiency
- Integration with quantum chemistry software (via interface)

## Usage Examples

### C++ Example

```cpp
#include "libfrag/libfrag.hpp"

// Configure MBE calculation
auto config = libfrag::MBEConfig::default_2body();
config.set_qm_method("B3LYP");
config.set_basis_set("def2-TZVP");

// Create molecule and calculator
libfrag::Molecule molecule(atoms);
auto qm_calculator = std::make_unique<MyQMCalculator>();
libfrag::MBECalculator calculator(std::move(qm_calculator), config);

// Perform calculation
auto results = calculator.calculate(molecule);
std::cout << results.summary_report() << std::endl;
```

### Python Example

```python
import libfrag
from libfrag.mbe import *

# Configure calculation
config = MBEConfig.default_2body()
config.set_qm_method("B3LYP")
config.set_basis_set("def2-TZVP")

# Create calculator and run
molecule = libfrag.Molecule(atoms)
calculator = MBECalculator(qm_interface, config)
results = calculator.calculate(molecule)

print(results.summary_report())
```

## Fragmentation Schemes

### 1. Atomic Fragmentation
Each atom becomes a separate fragment. Useful for:
- Very detailed analysis
- Small systems
- Understanding atomic contributions

### 2. Molecular Fragmentation
Connected components become fragments. Good for:
- Intermolecular interaction studies
- Molecular clusters
- Default choice for many systems

### 3. Functional Group Fragmentation
Chemical functional groups become fragments. Ideal for:
- Understanding chemical reactivity
- Large organic molecules
- Structure-activity relationships

### 4. Distance-Based Fragmentation
Spatial clustering creates fragments. Useful for:
- Large systems
- Adaptive fragmentation
- Reducing computational cost

### 5. Custom Fragmentation
User-defined fragment definitions. Allows:
- Expert knowledge incorporation
- Specialized applications
- Fine-grained control

## Integration with QM Software

The MBE module uses an abstract `QMCalculatorInterface` that can be implemented for different quantum chemistry packages:

- **PySCF**: Python-based quantum chemistry
- **Gaussian**: Commercial quantum chemistry suite
- **ORCA**: Free quantum chemistry program
- **PSI4**: Open-source quantum chemistry package
- **Q-Chem**: Commercial quantum chemistry software

Example implementation structure:

```cpp
class PySCFCalculator : public QMCalculatorInterface {
public:
    FragmentCalculationResult calculate_energy(
        const Fragment& fragment,
        const std::string& method,
        const std::string& basis_set,
        int charge, int multiplicity,
        const std::unordered_map<std::string, std::string>& options) override {
        
        // Convert fragment to PySCF input
        // Run PySCF calculation
        // Parse results and return
    }
};
```

## Performance Considerations

### Computational Scaling
The number of fragment combinations scales combinatorially:
- N fragments, 2-body: C(N,2) = N(N-1)/2 calculations
- N fragments, 3-body: C(N,3) = N(N-1)(N-2)/6 calculations

### Optimization Strategies
1. **Distance Filtering**: Skip distant fragment pairs
2. **Energy Thresholding**: Neglect small contributions
3. **Parallel Execution**: Run calculations in parallel
4. **Caching**: Avoid duplicate calculations
5. **Incremental Calculations**: Add higher-order terms as needed

### Memory Management
- Fragment results can be large (density matrices, MO coefficients)
- Streaming results to disk for very large calculations
- Configurable caching policies

## Error Analysis and Convergence

### Truncation Error Estimation
- Based on magnitude of highest-order terms
- Extrapolation methods for series convergence
- Statistical analysis of energy contributions

### Convergence Criteria
- Energy convergence thresholds
- Relative contribution limits
- Distance-based cutoffs

## Future Enhancements

### Planned Features
1. **Advanced Fragmentation**: Machine learning-based fragment selection
2. **Embedding Methods**: QM/MM and charge embedding
3. **Property Calculations**: Beyond energy (gradients, Hessians, properties)
4. **Visualization**: 3D visualization of fragment interactions
5. **Benchmark Database**: Standard test cases and reference results

### Integration Opportunities
- **Workflow Management**: Integration with scientific workflow systems
- **Cloud Computing**: Distributed calculation support
- **Database Storage**: Results database with search capabilities
- **Machine Learning**: ML-accelerated energy predictions

## Implementation Notes

This is scaffolding code that provides the structure and interfaces for MBE functionality. Key areas requiring implementation:

1. **Fragment Creation**: Converting atom indices to Fragment objects
2. **QM Integration**: Actual quantum chemistry software connections
3. **Functional Groups**: Pattern recognition for chemical groups
4. **Optimization**: Performance tuning and memory management
5. **Validation**: Comprehensive testing against reference calculations

The design emphasizes:
- **Modularity**: Clean separation of concerns
- **Extensibility**: Easy addition of new features
- **Performance**: Consideration for large-scale calculations
- **Usability**: Both expert and novice-friendly interfaces
- **Standards**: Compatibility with quantum chemistry conventions
