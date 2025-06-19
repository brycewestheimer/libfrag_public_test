# LibFrag Atom Class

The `Atom` class is a comprehensive representation of atoms designed specifically for quantum chemistry applications and seamless integration with MolSSI tools like QCEngine and QCElemental.

## Features

### Core Properties
- **Atomic identification**: Atomic number, element symbol, mass
- **Coordinates**: 3D position in Angstroms
- **Electronic properties**: Formal charge, partial charge, multiplicity
- **Quantum properties**: Extensible property storage system
- **Isotope support**: Mass number specification

### Key Capabilities
- **Multiple construction methods**: From atomic numbers, symbols, coordinate arrays
- **Geometric calculations**: Distance and vector calculations between atoms
- **Format compatibility**: Direct export to numpy arrays and QCElemental formats
- **Property management**: Store and retrieve quantum mechanical properties
- **Unit conversion**: Built-in constants for common QC unit conversions

## Quick Start

### C++ Usage

```cpp
#include "libfrag/atom.hpp"

// Create atoms
libfrag::Atom carbon(6, 0.0, 0.0, 0.0);
libfrag::Atom hydrogen("H", 1.089, 0.0, 0.0);

// Access properties
std::cout << "Symbol: " << carbon.symbol() << std::endl;
std::cout << "Mass: " << carbon.atomic_mass() << " amu" << std::endl;

// Calculate distance
double distance = carbon.distance_to(hydrogen);

// Set quantum properties
carbon.set_property("mulliken_charge", -0.123);
carbon.set_formal_charge(0.0);
```

### Python Usage

```python
import _libfrag as libfrag
import numpy as np

# Create atoms multiple ways
carbon = libfrag.Atom(6, 0.0, 0.0, 0.0)
nitrogen = libfrag.Atom("N", 1.5, 0.0, 0.0)
oxygen = libfrag.Atom(8, np.array([0.0, 1.5, 0.0]))

# Access properties
print(f"Carbon symbol: {carbon.symbol}")
print(f"Carbon mass: {carbon.atomic_mass} amu")

# Quantum properties
carbon.set_property("esp_charge", -0.145)
carbon.formal_charge = 0.0
carbon.partial_charge = -0.123

# Geometric calculations
distance = carbon.distance_to(nitrogen)
vector = carbon.vector_to(nitrogen)

# Export for QC software
coords_array = carbon.coordinates_array()  # numpy array
coords_list = carbon.coordinates_vector()  # Python list
```

## Integration with Quantum Chemistry Software

### QCElemental Compatibility

```python
# Create molecule data for QCElemental
atoms = [
    libfrag.Atom("C", 0.0, 0.0, 0.0),
    libfrag.Atom("H", 1.089, 0.0, 0.0),
    # ... more atoms
]

# Extract data in QCElemental format
geometry = np.array([atom.coordinates_vector() for atom in atoms]).flatten()
symbols = [atom.symbol for atom in atoms]
atomic_numbers = [atom.atomic_number for atom in atoms]

# Use with QCElemental
import qcelemental as qcel
molecule = qcel.models.Molecule(
    symbols=symbols,
    geometry=geometry,
    molecular_charge=0,
    molecular_multiplicity=1
)
```

### QCEngine Integration

The atom coordinates and properties can be directly used in QCEngine calculations:

```python
# Prepare input for QCEngine
input_data = {
    "molecule": {
        "symbols": [atom.symbol for atom in atoms],
        "geometry": [atom.coordinates_vector() for atom in atoms],
        "molecular_charge": sum(atom.formal_charge for atom in atoms),
        "molecular_multiplicity": 1
    },
    "driver": "energy",
    "model": {"method": "B3LYP", "basis": "6-31G*"}
}
```

## Quantum Properties

The Atom class supports arbitrary quantum properties through a key-value system:

```python
atom = libfrag.Atom("C", 0.0, 0.0, 0.0)

# Population analysis results
atom.set_property("mulliken_charge", -0.123)
atom.set_property("lowdin_charge", -0.098)
atom.set_property("esp_charge", -0.145)
atom.set_property("npa_charge", -0.134)

# Spin properties
atom.set_property("spin_density", 0.0)
atom.set_property("s_character", 0.25)
atom.set_property("p_character", 0.75)

# Energy contributions
atom.set_property("kinetic_energy", 37.123)
atom.set_property("potential_energy", -74.456)

# Access properties
charge = atom.get_property("mulliken_charge")  # Returns optional<double>
all_props = atom.properties  # Returns dictionary
```

## Unit Conversions

Built-in constants for common quantum chemistry unit conversions:

```python
# Available constants
libfrag.BOHR_TO_ANGSTROM  # 0.52917721067
libfrag.ANGSTROM_TO_BOHR  # 1.8897261246
libfrag.HARTREE_TO_EV     # 27.211386245988
libfrag.EV_TO_HARTREE     # 0.036749322176

# Example usage
coords_bohr = np.array(atom.coordinates_vector()) * libfrag.ANGSTROM_TO_BOHR
```

## Element Data

Comprehensive element data up to atomic number 118:

```python
# Static utility methods
z = libfrag.Atom.symbol_to_atomic_number("C")  # 6
symbol = libfrag.Atom.atomic_number_to_symbol(6)  # "C"

# Validation
is_valid = libfrag.Atom.is_valid_element_symbol("Xx")  # False
is_valid_z = libfrag.Atom.is_valid_atomic_number(200)  # False

# Automatic mass assignment
carbon = libfrag.Atom("C", 0.0, 0.0, 0.0)
print(carbon.atomic_mass)  # 12.011 amu
print(carbon.mass_number)  # 12 (most common isotope)
```

## Building and Installation

### Prerequisites

- C++17 compatible compiler
- CMake 3.15+
- xtensor and xtensor-python
- pybind11
- Python 3.6+ with NumPy

### Build Instructions

```bash
# Configure
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build --parallel

# Install (optional)
cmake --install build
```

### Python Package

```bash
# Install from source
pip install .

# Development install
pip install -e .
```

## API Reference

### Constructors

- `Atom()`  - Default constructor (hydrogen at origin)
- `Atom(int atomic_number, double x, double y, double z)` - From atomic number and coordinates
- `Atom(const std::string& symbol, double x, double y, double z)` - From symbol and coordinates
- `Atom(int atomic_number, const std::array<double, 3>& coords)` - From atomic number and array
- `Atom(int atomic_number, const xt::xarray<double>& coords)` - From numpy array (Python)

### Properties (Read-only)

- `atomic_number()` - Atomic number (Z)
- `symbol()` - Element symbol
- `coordinates()` - 3D coordinates array
- `x()`, `y()`, `z()` - Individual coordinates
- `atomic_mass()` - Atomic mass in amu
- `mass_number()` - Mass number of isotope

### Properties (Read-write)

- `formal_charge` - Formal charge
- `partial_charge` - Partial charge
- `multiplicity` - Spin multiplicity

### Methods

- `distance_to(other)` - Distance to another atom
- `vector_to(other)` - Vector to another atom
- `set_coordinates(...)` - Update coordinates
- `get_property(key)` - Get quantum property
- `set_property(key, value)` - Set quantum property
- `coordinates_array()` - Export as numpy array
- `coordinates_vector()` - Export as vector/list
- `to_string()` - String representation
- `to_xyz_string()` - XYZ format string

## Examples

See the included examples:
- `python/tests/test_atom.py` - Comprehensive test suite
- `python/notebooks/atom_demo.ipynb` - Interactive Jupyter notebook
- `examples/` - C++ usage examples

## License

This project is licensed under the terms specified in the LICENCE.txt file.
