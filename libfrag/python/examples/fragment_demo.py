#!/usr/bin/env python3
"""
Simple example demonstrating the Fragment class Python API
"""

import sys
import os

# Add the build directory to the path for importing the module
# This would typically be handled by proper installation
try:
    import libfrag
    print("✓ Successfully imported libfrag module")
except ImportError as e:
    print(f"✗ Failed to import libfrag module: {e}")
    print("Make sure the module is built and in your Python path")
    sys.exit(1)

def main():
    print("=" * 60)
    print("Fragment Class Python API Demo")
    print("=" * 60)
    
    # Create some atoms
    print("\n1. Creating atoms...")
    atoms = [
        libfrag.Atom(6, 0.0, 0.0, 0.0),    # Carbon
        libfrag.Atom(1, 1.1, 0.0, 0.0),    # Hydrogen
        libfrag.Atom(1, -0.5, 0.9, 0.0),   # Hydrogen
        libfrag.Atom(1, -0.5, -0.9, 0.0),  # Hydrogen
        libfrag.Atom(8, 5.0, 0.0, 0.0),    # Oxygen (separate)
        libfrag.Atom(1, 5.5, 0.8, 0.0),    # Hydrogen
    ]
    
    print(f"Created {len(atoms)} atoms:")
    for i, atom in enumerate(atoms):
        print(f"  {i}: {atom.symbol} at ({atom.x:.1f}, {atom.y:.1f}, {atom.z:.1f})")
    
    # Create bonds (methane + water, disconnected)
    print("\n2. Creating bonds...")
    bonds = [(0, 1), (0, 2), (0, 3), (4, 5)]  # Methane + water bonds
    print(f"Bond connectivity: {bonds}")
    
    # Create fragment
    print("\n3. Creating fragment...")
    fragment = libfrag.Fragment(atoms, bonds, "complex_molecule")
    
    print(f"Fragment ID: {fragment.fragment_id}")
    print(f"Atoms: {fragment.atom_count}")
    print(f"Bonds: {fragment.bond_count}")
    print(f"Is connected: {fragment.is_connected()}")
    print(f"Generation level: {fragment.generation_level}")
    print(f"Is root fragment: {fragment.is_root_fragment()}")
    
    # Calculate properties
    print("\n4. Fragment properties...")
    properties = fragment.calculate_fragment_properties()
    for name, value in properties.items():
        print(f"  {name}: {value}")
    
    # Get composition
    print("\n5. Fragment composition...")
    composition = fragment.fragment_composition()
    for element, count in composition.items():
        print(f"  {element}: {count}")
    
    # Fragment by connectivity
    print("\n6. Fragmenting by connectivity...")
    subfragments = fragment.fragment_by_connectivity()
    print(f"Created {len(subfragments)} subfragments:")
    
    for i, subfrag in enumerate(subfragments):
        comp = subfrag.fragment_composition()
        comp_str = ", ".join([f"{elem}{count}" for elem, count in comp.items()])
        print(f"  Subfragment {i}: ID='{subfrag.fragment_id}', atoms={subfrag.atom_count}, composition={comp_str}")
    
    # Display fragment tree
    print("\n7. Fragment tree structure...")
    print(fragment.to_tree_string())
    
    # Test fragment links
    if len(subfragments) >= 2:
        print("\n8. Testing fragment links...")
        
        # Create a link between first two subfragments
        link = libfrag.FragmentLink(
            subfragments[0].fragment_id,
            subfragments[1].fragment_id,
            0, 0,  # Connect first atoms
            libfrag.BondType.SINGLE,
            1.0
        )
        
        subfragments[0].add_fragment_link(link)
        print(f"Added link: {link}")
        
        linked_ids = subfragments[0].linked_fragment_ids
        print(f"Fragment '{subfragments[0].fragment_id}' is linked to: {list(linked_ids)}")
    
    # Test factory methods
    print("\n9. Testing factory methods...")
    
    # Create molecule first
    molecule = libfrag.Molecule(atoms[:4])  # Just the methane part
    factory_fragment = libfrag.Fragment.create_from_molecule(molecule, "factory_test")
    print(f"Factory fragment: {factory_fragment.fragment_id}, atoms: {factory_fragment.atom_count}")
    
    # Test hierarchical creation
    hierarchy_fragment = libfrag.Fragment.create_fragment_hierarchy(
        molecule, max_fragment_size=2, id_prefix="hierarchy"
    )
    print(f"Hierarchy fragment: {hierarchy_fragment.fragment_id}")
    print(f"  Subfragments: {hierarchy_fragment.subfragment_count}")
    
    # Test utility functions
    print("\n10. Testing utility functions...")
    
    test_fragment = libfrag.create_test_fragment()
    print(f"Test fragment: {test_fragment.fragment_id}, atoms: {test_fragment.atom_count}")
    
    complex_test = libfrag.create_complex_test_fragment()
    print(f"Complex test fragment: {complex_test.fragment_id}, atoms: {complex_test.atom_count}")
    
    # Test fragmentation on complex test
    complex_subfragments = complex_test.fragment_by_connectivity()
    print(f"Complex test fragmented into {len(complex_subfragments)} pieces")
    
    print("\n" + "=" * 60)
    print("Fragment Class Demo Complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
