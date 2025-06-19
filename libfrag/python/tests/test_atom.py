"""
Test and demonstration of the libfrag Atom class.

This module shows how to use the Atom class for quantum chemistry applications,
with compatibility for MolSSI tools like QCEngine and QCElemental.
"""

import numpy as np
import sys
import os

# Add the build directory to the path (adjust as needed)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'python', 'src'))

try:
    import _libfrag as libfrag
    print(f"Successfully imported libfrag version {libfrag.__version__}")
except ImportError as e:
    print(f"Failed to import libfrag: {e}")
    print("Make sure the library has been built and is in the Python path")
    sys.exit(1)


def test_basic_atom_creation():
    """Test basic atom creation methods."""
    print("\n=== Basic Atom Creation ===")
    
    # Create atoms different ways
    carbon = libfrag.Atom(6, 0.0, 0.0, 0.0)
    nitrogen = libfrag.Atom("N", 1.5, 0.0, 0.0)
    oxygen = libfrag.Atom(8, np.array([0.0, 1.5, 0.0]))
    
    print(f"Carbon: {carbon}")
    print(f"Nitrogen: {nitrogen}")
    print(f"Oxygen: {oxygen}")
    
    # Test properties
    print(f"\nCarbon atomic number: {carbon.atomic_number}")
    print(f"Carbon symbol: {carbon.symbol}")
    print(f"Carbon atomic mass: {carbon.atomic_mass} amu")
    print(f"Carbon coordinates: {carbon.coordinates}")


def test_quantum_properties():
    """Test quantum property management."""
    print("\n=== Quantum Properties ===")
    
    atom = libfrag.Atom("C", 0.0, 0.0, 0.0)
    
    # Set various quantum properties
    atom.set_property("mulliken_charge", -0.123)
    atom.set_property("esp_charge", -0.145)
    atom.set_property("spin_density", 0.0)
    atom.formal_charge = 0.0
    atom.partial_charge = -0.123
    atom.multiplicity = 1
    
    print(f"Mulliken charge: {atom.get_property('mulliken_charge')}")
    print(f"ESP charge: {atom.get_property('esp_charge')}")
    print(f"Formal charge: {atom.formal_charge}")
    print(f"Partial charge: {atom.partial_charge}")
    print(f"All properties: {atom.properties}")


def test_molecular_geometry():
    """Test geometric calculations between atoms."""
    print("\n=== Molecular Geometry ===")
    
    # Create a simple water molecule
    oxygen = libfrag.Atom("O", 0.0, 0.0, 0.0)
    hydrogen1 = libfrag.Atom("H", 0.757, 0.586, 0.0)
    hydrogen2 = libfrag.Atom("H", -0.757, 0.586, 0.0)
    
    print("Water molecule coordinates:")
    print(f"O: {oxygen.to_xyz_string()}")
    print(f"H: {hydrogen1.to_xyz_string()}")
    print(f"H: {hydrogen2.to_xyz_string()}")
    
    # Calculate distances
    oh1_distance = oxygen.distance_to(hydrogen1)
    oh2_distance = oxygen.distance_to(hydrogen2)
    hh_distance = hydrogen1.distance_to(hydrogen2)
    
    print(f"\nBond distances (Angstroms):")
    print(f"O-H1: {oh1_distance:.3f}")
    print(f"O-H2: {oh2_distance:.3f}")
    print(f"H1-H2: {hh_distance:.3f}")
    
    # Calculate vectors
    oh1_vector = oxygen.vector_to(hydrogen1)
    print(f"\nO -> H1 vector: {oh1_vector}")


def test_qc_compatibility():
    """Test features useful for quantum chemistry software integration."""
    print("\n=== QC Software Compatibility ===")
    
    # Unit conversion constants
    print(f"Bohr to Angstrom: {libfrag.BOHR_TO_ANGSTROM}")
    print(f"Angstrom to Bohr: {libfrag.ANGSTROM_TO_BOHR}")
    print(f"Hartree to eV: {libfrag.HARTREE_TO_EV}")
    
    # Create atoms for a molecule
    atoms = [
        libfrag.Atom("C", 0.0, 0.0, 0.0),
        libfrag.Atom("H", 1.089, 0.0, 0.0),
        libfrag.Atom("H", -0.363, 1.027, 0.0),
        libfrag.Atom("H", -0.363, -0.513, 0.889),
        libfrag.Atom("H", -0.363, -0.513, -0.889)
    ]
    
    print("\nMethane molecule:")
    for i, atom in enumerate(atoms):
        print(f"{i}: {atom.to_xyz_string()}")
    
    # Convert to formats useful for QC packages
    print("\nCoordinates as numpy arrays:")
    for atom in atoms:
        coords_array = atom.coordinates_array()
        print(f"{atom.symbol}: {coords_array}")
    
    # Create coordinate matrix (useful for QCElemental)
    geometry_matrix = np.array([atom.coordinates_vector() for atom in atoms])
    symbols = [atom.symbol for atom in atoms]
    atomic_numbers = [atom.atomic_number for atom in atoms]
    
    print(f"\nGeometry matrix shape: {geometry_matrix.shape}")
    print(f"Symbols: {symbols}")
    print(f"Atomic numbers: {atomic_numbers}")


def test_static_methods():
    """Test static utility methods."""
    print("\n=== Static Utility Methods ===")
    
    # Symbol/atomic number conversion
    print(f"'C' -> {libfrag.Atom.symbol_to_atomic_number('C')}")
    print(f"6 -> '{libfrag.Atom.atomic_number_to_symbol(6)}'")
    
    # Validation
    print(f"Is 'C' valid? {libfrag.Atom.is_valid_element_symbol('C')}")
    print(f"Is 'Xx' valid? {libfrag.Atom.is_valid_element_symbol('Xx')}")
    print(f"Is Z=6 valid? {libfrag.Atom.is_valid_atomic_number(6)}")
    print(f"Is Z=200 valid? {libfrag.Atom.is_valid_atomic_number(200)}")


def create_qcelemental_compatible_molecule():
    """Create a molecule data structure compatible with QCElemental format."""
    print("\n=== QCElemental Compatible Structure ===")
    
    # Create benzene molecule
    benzene_coords = [
        [0.0000,  1.3970, 0.0000],  # C
        [1.2098,  0.6985, 0.0000],  # C
        [1.2098, -0.6985, 0.0000],  # C
        [0.0000, -1.3970, 0.0000],  # C
        [-1.2098, -0.6985, 0.0000], # C
        [-1.2098,  0.6985, 0.0000], # C
        [0.0000,  2.4820, 0.0000],  # H
        [2.1498,  1.2410, 0.0000],  # H
        [2.1498, -1.2410, 0.0000],  # H
        [0.0000, -2.4820, 0.0000],  # H
        [-2.1498, -1.2410, 0.0000], # H
        [-2.1498,  1.2410, 0.0000]  # H
    ]
    
    benzene_symbols = ['C'] * 6 + ['H'] * 6
    
    # Create libfrag atoms
    atoms = []
    for symbol, coords in zip(benzene_symbols, benzene_coords):
        atom = libfrag.Atom(symbol, coords[0], coords[1], coords[2])
        atoms.append(atom)
    
    # Extract data in QCElemental format
    geometry = np.array([atom.coordinates_vector() for atom in atoms]).flatten()
    symbols = [atom.symbol for atom in atoms]
    atomic_numbers = [atom.atomic_number for atom in atoms]
    masses = [atom.atomic_mass for atom in atoms]
    
    print("Benzene molecule (QCElemental compatible):")
    print(f"Symbols: {symbols}")
    print(f"Atomic numbers: {atomic_numbers}")
    print(f"Geometry (flattened): shape={geometry.shape}")
    print(f"Masses: {masses}")
    
    # This could be used directly with QCElemental Molecule class:
    # mol = qcel.models.Molecule(
    #     symbols=symbols,
    #     geometry=geometry,
    #     molecular_charge=0,
    #     molecular_multiplicity=1
    # )
    
    return {
        'symbols': symbols,
        'geometry': geometry,
        'atomic_numbers': atomic_numbers,
        'masses': masses,
        'molecular_charge': 0,
        'molecular_multiplicity': 1
    }


if __name__ == "__main__":
    print("LibFrag Atom Class Test Suite")
    print("=" * 50)
    
    try:
        test_basic_atom_creation()
        test_quantum_properties()
        test_molecular_geometry()
        test_qc_compatibility()
        test_static_methods()
        mol_data = create_qcelemental_compatible_molecule()
        
        print("\n" + "=" * 50)
        print("All tests completed successfully!")
        print("\nThe Atom class is ready for use with:")
        print("- QCElemental for molecule data structures")
        print("- QCEngine for quantum chemistry calculations")
        print("- MolSSI Driver Interface for workflow management")
        
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
