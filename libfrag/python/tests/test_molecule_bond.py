"""
Test and demonstration of the libfrag Bond and Molecule classes.

This module shows how to use the Bond and Molecule classes for quantum chemistry
applications, with compatibility for MolSSI tools like QCEngine and QCElemental.
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


def test_bond_creation():
    """Test basic bond creation and properties."""
    print("\n=== Bond Creation and Properties ===")
    
    # Create atoms
    carbon = libfrag.Atom(6, 0.0, 0.0, 0.0)
    hydrogen = libfrag.Atom("H", 1.089, 0.0, 0.0)
    nitrogen = libfrag.Atom("N", 0.0, 1.5, 0.0)
    
    # Create different types of bonds
    ch_bond = libfrag.Bond(carbon, hydrogen, libfrag.BondType.SINGLE)
    cn_bond = libfrag.Bond(carbon, nitrogen, 1.5)  # Partial double bond
    
    print(f"C-H bond: {ch_bond}")
    print(f"C-N bond: {cn_bond}")
    
    # Test bond properties
    print(f"\nC-H bond length: {ch_bond.bond_length:.3f} Angstroms")
    print(f"C-H bond order: {ch_bond.bond_order}")
    print(f"C-H is polar: {ch_bond.is_polar_bond()}")
    print(f"C-H is strong bond: {ch_bond.is_strong_bond()}")
    
    # Set quantum properties
    ch_bond.set_property("wiberg_bond_order", 0.98)
    ch_bond.set_property("mayer_bond_order", 0.95)
    ch_bond.bond_strength = 413.0  # C-H bond strength in kcal/mol
    
    print(f"\nC-H Wiberg bond order: {ch_bond.get_property('wiberg_bond_order')}")
    print(f"C-H bond strength: {ch_bond.bond_strength} kcal/mol")
    print(f"All C-H properties: {ch_bond.properties}")


def test_molecule_creation():
    """Test basic molecule creation methods."""
    print("\n=== Molecule Creation ===")
    
    # Create empty molecule and add atoms
    mol = libfrag.Molecule()
    print(f"Empty molecule: {mol}")
    
    # Add atoms individually
    c_idx = mol.add_atom(6, 0.0, 0.0, 0.0)
    h1_idx = mol.add_atom("H", 1.089, 0.0, 0.0)
    h2_idx = mol.add_atom("H", -0.363, 1.027, 0.0)
    h3_idx = mol.add_atom("H", -0.363, -0.513, 0.889)
    h4_idx = mol.add_atom("H", -0.363, -0.513, -0.889)
    
    print(f"Molecule with atoms: {mol}")
    print(f"Atom count: {mol.atom_count}")
    
    # Add bonds
    mol.add_bond(c_idx, h1_idx, libfrag.BondType.SINGLE)
    mol.add_bond(c_idx, h2_idx, libfrag.BondType.SINGLE)
    mol.add_bond(c_idx, h3_idx, libfrag.BondType.SINGLE)
    mol.add_bond(c_idx, h4_idx, libfrag.BondType.SINGLE)
    
    print(f"Methane molecule: {mol}")
    print(f"Bond count: {mol.bond_count}")


def test_molecule_from_atoms():
    """Test creating molecule from list of atoms."""
    print("\n=== Molecule from Atom List ===")
    
    # Create water molecule
    atoms = [
        libfrag.Atom("O", 0.0, 0.0, 0.0),
        libfrag.Atom("H", 0.757, 0.586, 0.0),
        libfrag.Atom("H", -0.757, 0.586, 0.0)
    ]
    
    # Define connectivity (atom index pairs)
    bonds = [(0, 1), (0, 2)]  # O-H bonds
    bond_types = [libfrag.BondType.SINGLE, libfrag.BondType.SINGLE]
    
    water = libfrag.Molecule(atoms, bonds, bond_types)
    print(f"Water molecule: {water}")
    print(f"Water is connected: {water.is_connected()}")
    
    # Test molecular properties
    print(f"\nWater molecular mass: {water.molecular_mass():.3f} amu")
    print(f"Water molecular charge: {water.molecular_charge()}")
    print(f"Water center of mass: {water.center_of_mass()}")
    
    # Test coordination numbers
    print(f"O coordination number: {water.coordination_number(0)}")
    print(f"H coordination number: {water.coordination_number(1)}")
    
    # Test bonded atoms
    bonded_to_oxygen = water.bonded_atoms(0)
    print(f"Atoms bonded to oxygen: {bonded_to_oxygen}")


def test_molecular_geometry():
    """Test molecular geometry operations."""
    print("\n=== Molecular Geometry Operations ===")
    
    # Create benzene molecule (simplified)
    benzene_atoms = [
        libfrag.Atom("C", 0.0000,  1.3970, 0.0000),
        libfrag.Atom("C", 1.2098,  0.6985, 0.0000),
        libfrag.Atom("C", 1.2098, -0.6985, 0.0000),
        libfrag.Atom("C", 0.0000, -1.3970, 0.0000),
        libfrag.Atom("C", -1.2098, -0.6985, 0.0000),
        libfrag.Atom("C", -1.2098,  0.6985, 0.0000)
    ]
    
    # Define aromatic ring connectivity
    benzene_bonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
    benzene_bond_types = [libfrag.BondType.AROMATIC] * 6
    
    benzene = libfrag.Molecule(benzene_atoms, benzene_bonds, benzene_bond_types)
    print(f"Benzene: {benzene}")
    
    # Test geometry properties
    print(f"Benzene diameter: {benzene.molecular_diameter():.3f} Angstroms")
    print(f"Benzene bounding box: {benzene.bounding_box()}")
    
    # Center at origin
    original_com = benzene.center_of_mass()
    print(f"Original center of mass: {original_com}")
    
    benzene.center_at_origin()
    new_com = benzene.center_of_mass()
    print(f"New center of mass: {new_com}")
    
    # Translate molecule
    benzene.translate([1.0, 2.0, 3.0])
    translated_com = benzene.center_of_mass()
    print(f"Translated center of mass: {translated_com}")


def test_bond_analysis():
    """Test bond analysis and properties."""
    print("\n=== Bond Analysis ===")
    
    # Create ethane molecule
    atoms = [
        libfrag.Atom("C", 0.0, 0.0, 0.0),
        libfrag.Atom("C", 1.54, 0.0, 0.0),
        libfrag.Atom("H", -0.51, 1.02, 0.0),
        libfrag.Atom("H", -0.51, -0.51, 0.88),
        libfrag.Atom("H", -0.51, -0.51, -0.88),
        libfrag.Atom("H", 2.05, 1.02, 0.0),
        libfrag.Atom("H", 2.05, -0.51, 0.88),
        libfrag.Atom("H", 2.05, -0.51, -0.88)
    ]
    
    ethane = libfrag.Molecule(atoms)
    
    # Add C-C bond
    cc_bond_idx = ethane.add_bond(0, 1, libfrag.BondType.SINGLE)
    
    # Add C-H bonds
    for h_idx in range(2, 5):
        ethane.add_bond(0, h_idx, libfrag.BondType.SINGLE)
    for h_idx in range(5, 8):
        ethane.add_bond(1, h_idx, libfrag.BondType.SINGLE)
    
    print(f"Ethane: {ethane}")
    
    # Analyze C-C bond
    cc_bond = ethane.bond(cc_bond_idx)
    print(f"\nC-C bond: {cc_bond}")
    print(f"C-C bond length: {cc_bond.bond_length:.3f} Angstroms")
    print(f"C-C involves carbon: {cc_bond.involves_atom_type(6)}")
    print(f"C-C involves hydrogen: {cc_bond.involves_atom_type(1)}")
    
    # Set bond properties
    cc_bond.set_property("rotation_barrier", 2.9)  # kcal/mol
    cc_bond.equilibrium_bond_length = 1.54
    cc_bond.bond_strength = 83.0  # kcal/mol
    
    print(f"C-C rotation barrier: {cc_bond.get_property('rotation_barrier')} kcal/mol")
    print(cc_bond.to_detailed_string())


def test_qc_compatibility():
    """Test quantum chemistry software compatibility."""
    print("\n=== QC Software Compatibility ===")
    
    # Create a molecule for QC calculations
    mol = libfrag.Molecule.create_linear_molecule(["C", "N", "O"], 1.2)
    print(f"Linear CNO molecule: {mol}")
    
    # Get data in QC format
    atomic_numbers = mol.atomic_numbers()
    symbols = mol.element_symbols()
    geometry = mol.geometry_array()
    connectivity = mol.connectivity_matrix()
    
    print(f"Atomic numbers: {atomic_numbers}")
    print(f"Element symbols: {symbols}")
    print(f"Geometry array shape: {geometry.shape}")
    print(f"Connectivity matrix shape: {connectivity.shape}")
    
    # Set molecular properties for QC calculations
    mol.molecular_multiplicity = 2  # Doublet
    mol.set_property("total_energy", -168.123456)  # Hartrees
    mol.set_property("homo_energy", -0.345)
    mol.set_property("lumo_energy", 0.123)
    
    print(f"\nMolecular multiplicity: {mol.molecular_multiplicity}")
    print(f"Total energy: {mol.get_property('total_energy')} hartrees")
    print(f"HOMO energy: {mol.get_property('homo_energy')} hartrees")
    
    # Generate XYZ format
    xyz_string = mol.to_xyz_string("CNO linear molecule for QC calculation")
    print(f"\nXYZ format:\n{xyz_string}")


def test_molecular_fragments():
    """Test molecular fragmentation analysis."""
    print("\n=== Molecular Fragmentation ===")
    
    # Create two separate fragments
    fragment1_atoms = [
        libfrag.Atom("C", 0.0, 0.0, 0.0),
        libfrag.Atom("H", 1.0, 0.0, 0.0)
    ]
    
    fragment2_atoms = [
        libfrag.Atom("N", 5.0, 0.0, 0.0),  # Far away
        libfrag.Atom("H", 6.0, 0.0, 0.0)
    ]
    
    # Combine into one molecule
    all_atoms = fragment1_atoms + fragment2_atoms
    bonds = [(0, 1), (2, 3)]  # Only bonds within fragments
    
    mol = libfrag.Molecule(all_atoms, bonds)
    print(f"Disconnected molecule: {mol}")
    print(f"Is connected: {mol.is_connected()}")
    
    # Get separate fragments
    fragments = mol.get_fragments()
    print(f"Number of fragments: {len(fragments)}")
    
    for i, fragment in enumerate(fragments):
        print(f"Fragment {i}: {fragment}")
        print(f"  Atoms: {fragment.atom_count}")
        print(f"  Bonds: {fragment.bond_count}")


def test_xyz_file_io():
    """Test XYZ file format I/O."""
    print("\n=== XYZ File I/O ===")
    
    # Create sample XYZ string
    xyz_content = """3
Water molecule
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000"""
    
    # Parse from XYZ
    water = libfrag.Molecule.from_xyz_string(xyz_content)
    print(f"Water from XYZ: {water}")
    
    # Convert back to XYZ
    xyz_output = water.to_xyz_string("Reconstructed water molecule")
    print(f"XYZ output:\n{xyz_output}")


def create_qcelemental_molecule():
    """Create QCElemental-compatible molecule data."""
    print("\n=== QCElemental Compatible Structure ===")
    
    # Create caffeine molecule (simplified structure)
    caffeine_mol = libfrag.Molecule()
    
    # Add atoms (simplified caffeine structure)
    atom_data = [
        ("C", 0.000, 0.000, 0.000),
        ("N", 1.350, 0.000, 0.000),
        ("C", 2.000, 1.200, 0.000),
        ("C", 1.300, 2.400, 0.000),
        ("C", -0.100, 2.400, 0.000),
        ("N", -0.750, 1.200, 0.000),
        ("O", 3.220, 1.200, 0.000),
        ("N", 2.000, 3.600, 0.000)
    ]
    
    for symbol, x, y, z in atom_data:
        caffeine_mol.add_atom(libfrag.Atom(symbol, x, y, z))
    
    # Add some bonds (simplified)
    bonds_to_add = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0), (2,6), (3,7)]
    for i, j in bonds_to_add:
        caffeine_mol.add_bond(i, j)
    
    print(f"Caffeine-like molecule: {caffeine_mol}")
    
    # Extract QCElemental-compatible data
    qcel_data = {
        'symbols': caffeine_mol.element_symbols(),
        'geometry': caffeine_mol.geometry_array(),
        'atomic_numbers': caffeine_mol.atomic_numbers(),
        'molecular_charge': caffeine_mol.molecular_charge(),
        'molecular_multiplicity': caffeine_mol.molecular_multiplicity,
        'connectivity': caffeine_mol.connectivity_matrix()
    }
    
    print(f"QCElemental data prepared:")
    print(f"  Symbols: {qcel_data['symbols']}")
    print(f"  Geometry shape: {qcel_data['geometry'].shape}")
    print(f"  Connectivity shape: {qcel_data['connectivity'].shape}")
    
    return qcel_data


if __name__ == "__main__":
    print("LibFrag Bond and Molecule Classes Test Suite")
    print("=" * 60)
    
    try:
        test_bond_creation()
        test_molecule_creation()
        test_molecule_from_atoms()
        test_molecular_geometry()
        test_bond_analysis()
        test_qc_compatibility()
        test_molecular_fragments()
        test_xyz_file_io()
        qcel_data = create_qcelemental_molecule()
        
        print("\n" + "=" * 60)
        print("All tests completed successfully!")
        print("\nThe Bond and Molecule classes are ready for use with:")
        print("- QCElemental for molecule data structures")
        print("- QCEngine for quantum chemistry calculations")
        print("- MolSSI Driver Interface for workflow management")
        print("- Fragment-based quantum chemistry methods")
        
    except Exception as e:
        print(f"\nTest failed with error: {e}")
        import traceback
        traceback.print_exc()
