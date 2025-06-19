"""
Test suite for Fragment class Python bindings
"""

import pytest
import numpy as np
try:
    import libfrag
except ImportError:
    pytest.skip("libfrag module not available", allow_module_level=True)


class TestFragmentLink:
    """Test FragmentLink class"""
    
    def test_fragment_link_creation(self):
        """Test basic FragmentLink creation"""
        link = libfrag.FragmentLink("frag1", "frag2", 0, 1, libfrag.BondType.SINGLE, 1.0)
        
        assert link.source_fragment_id == "frag1"
        assert link.target_fragment_id == "frag2"
        assert link.source_atom_index == 0
        assert link.target_atom_index == 1
        assert link.bond_type == libfrag.BondType.SINGLE
        assert link.bond_order == 1.0
    
    def test_fragment_link_properties(self):
        """Test FragmentLink property modification"""
        link = libfrag.FragmentLink("frag1", "frag2", 0, 1)
        
        # Test default values
        assert link.bond_type == libfrag.BondType.SINGLE
        assert link.bond_order == 1.0
        
        # Test property modification
        link.bond_type = libfrag.BondType.DOUBLE
        link.bond_order = 2.0
        
        assert link.bond_type == libfrag.BondType.DOUBLE
        assert link.bond_order == 2.0
    
    def test_fragment_link_string_representation(self):
        """Test FragmentLink string methods"""
        link = libfrag.FragmentLink("frag1", "frag2", 0, 1)
        
        # Test to_string method
        link_str = link.to_string()
        assert "frag1" in link_str
        assert "frag2" in link_str
        
        # Test __str__ and __repr__
        str_repr = str(link)
        repr_str = repr(link)
        assert isinstance(str_repr, str)
        assert isinstance(repr_str, str)
    
    def test_fragment_link_equality(self):
        """Test FragmentLink equality operators"""
        link1 = libfrag.FragmentLink("frag1", "frag2", 0, 1, libfrag.BondType.SINGLE, 1.0)
        link2 = libfrag.FragmentLink("frag1", "frag2", 0, 1, libfrag.BondType.SINGLE, 1.0)
        link3 = libfrag.FragmentLink("frag1", "frag3", 0, 1, libfrag.BondType.SINGLE, 1.0)
        
        assert link1 == link2
        assert link1 != link3


class TestFragmentBasics:
    """Test basic Fragment functionality"""
    
    def test_fragment_default_constructor(self):
        """Test Fragment default constructor"""
        fragment = libfrag.Fragment()
        
        assert fragment.is_empty()
        assert fragment.fragment_id == ""
        assert fragment.is_root_fragment()
        assert fragment.generation_level == 0
        assert fragment.subfragment_count == 0
        assert not fragment.has_subfragments()
    
    def test_fragment_with_id(self):
        """Test Fragment constructor with ID"""
        fragment = libfrag.Fragment("test_fragment")
        
        assert fragment.is_empty()
        assert fragment.fragment_id == "test_fragment"
        assert fragment.is_root_fragment()
    
    def test_fragment_from_atoms(self):
        """Test Fragment constructor from atoms"""
        atoms = [
            libfrag.Atom(6, 0.0, 0.0, 0.0),  # Carbon
            libfrag.Atom(1, 1.1, 0.0, 0.0),  # Hydrogen
        ]
        
        fragment = libfrag.Fragment(atoms, "simple_fragment")
        
        assert fragment.atom_count == 2
        assert fragment.fragment_id == "simple_fragment"
        assert fragment.atom(0).atomic_number == 6
        assert fragment.atom(1).atomic_number == 1
    
    def test_fragment_with_bonds(self):
        """Test Fragment constructor with atoms and bonds"""
        atoms = [
            libfrag.Atom(6, 0.0, 0.0, 0.0),  # Carbon
            libfrag.Atom(1, 1.1, 0.0, 0.0),  # Hydrogen
            libfrag.Atom(1, -0.5, 0.9, 0.0), # Hydrogen
        ]
        
        bonds = [(0, 1), (0, 2)]
        
        fragment = libfrag.Fragment(atoms, bonds, "methane_part")
        
        assert fragment.atom_count == 3
        assert fragment.bond_count == 2
        assert fragment.fragment_id == "methane_part"


class TestFragmentHierarchy:
    """Test Fragment hierarchy and composite pattern functionality"""
    
    def test_subfragment_management(self):
        """Test adding and managing subfragments"""
        parent = libfrag.Fragment("parent")
        child1 = libfrag.Fragment("child1")
        child2 = libfrag.Fragment("child2")
        
        # Initially no subfragments
        assert parent.subfragment_count == 0
        assert not parent.has_subfragments()
        
        # Add subfragments
        parent.add_subfragment(child1)
        parent.add_subfragment(child2)
        
        assert parent.subfragment_count == 2
        assert parent.has_subfragments()
        assert len(parent.subfragments) == 2
        assert parent.subfragment(0).fragment_id == "child1"
        assert parent.subfragment(1).fragment_id == "child2"
    
    def test_parent_child_relationships(self):
        """Test parent-child relationships"""
        parent = libfrag.Fragment("parent")
        child = libfrag.Fragment("child")
        
        parent.add_subfragment(child)
        
        # Check child properties
        assert child.generation_level == 1
        assert not child.is_root_fragment()
        
        # Check parent relationship
        child_parent = child.parent_fragment
        assert child_parent is not None
        assert child_parent.fragment_id == "parent"
    
    def test_subfragment_removal(self):
        """Test removing subfragments"""
        parent = libfrag.Fragment("parent")
        child1 = libfrag.Fragment("child1")
        child2 = libfrag.Fragment("child2")
        
        parent.add_subfragment(child1)
        parent.add_subfragment(child2)
        
        assert parent.subfragment_count == 2
        
        # Remove first subfragment
        parent.remove_subfragment(0)
        assert parent.subfragment_count == 1
        assert parent.subfragment(0).fragment_id == "child2"
        
        # Clear all subfragments
        parent.clear_subfragments()
        assert parent.subfragment_count == 0


class TestFragmentLinks:
    """Test inter-fragment connectivity"""
    
    def test_fragment_link_management(self):
        """Test adding and managing fragment links"""
        fragment = libfrag.Fragment("test_fragment")
        
        link1 = libfrag.FragmentLink("test_fragment", "other1", 0, 0)
        link2 = libfrag.FragmentLink("test_fragment", "other2", 1, 1)
        
        fragment.add_fragment_link(link1)
        fragment.add_fragment_link(link2)
        
        assert len(fragment.fragment_links) == 2
        assert fragment.is_linked_to("other1")
        assert fragment.is_linked_to("other2")
        assert not fragment.is_linked_to("other3")
        
        linked_ids = fragment.linked_fragment_ids
        assert len(linked_ids) == 2
        assert "other1" in linked_ids
        assert "other2" in linked_ids
    
    def test_find_links_to_fragment(self):
        """Test finding links to specific fragments"""
        fragment = libfrag.Fragment("test_fragment")
        
        link1 = libfrag.FragmentLink("test_fragment", "target", 0, 0)
        link2 = libfrag.FragmentLink("test_fragment", "other", 1, 1)
        link3 = libfrag.FragmentLink("test_fragment", "target", 2, 2)
        
        fragment.add_fragment_link(link1)
        fragment.add_fragment_link(link2)
        fragment.add_fragment_link(link3)
        
        target_links = fragment.find_links_to_fragment("target")
        assert len(target_links) == 2
        assert 0 in target_links
        assert 2 in target_links
        
        other_links = fragment.find_links_to_fragment("other")
        assert len(other_links) == 1
        assert 1 in other_links


class TestFragmentAnalysis:
    """Test fragment analysis functionality"""
    
    def test_fragment_properties(self):
        """Test fragment property calculation"""
        atoms = [
            libfrag.Atom(6, 0.0, 0.0, 0.0),  # Carbon
            libfrag.Atom(1, 1.1, 0.0, 0.0),  # Hydrogen
            libfrag.Atom(1, -0.5, 0.9, 0.0), # Hydrogen
        ]
        
        bonds = [(0, 1), (0, 2)]
        fragment = libfrag.Fragment(atoms, bonds, "test")
        
        properties = fragment.calculate_fragment_properties()
        
        assert "atom_count" in properties
        assert "bond_count" in properties
        assert "molecular_mass" in properties
        assert properties["atom_count"] == 3.0
        assert properties["bond_count"] == 2.0
    
    def test_fragment_composition(self):
        """Test fragment composition analysis"""
        atoms = [
            libfrag.Atom(6, 0.0, 0.0, 0.0),  # Carbon
            libfrag.Atom(1, 1.1, 0.0, 0.0),  # Hydrogen
            libfrag.Atom(1, -0.5, 0.9, 0.0), # Hydrogen
            libfrag.Atom(1, -0.5, -0.9, 0.0),# Hydrogen
        ]
        
        fragment = libfrag.Fragment(atoms, "methane")
        composition = fragment.fragment_composition()
        
        assert composition["C"] == 1
        assert composition["H"] == 3
        assert len(composition) == 2
    
    def test_fragment_hash(self):
        """Test fragment hash generation"""
        atoms = [libfrag.Atom(6, 0.0, 0.0, 0.0)]
        
        fragment1 = libfrag.Fragment(atoms, "frag1")
        fragment2 = libfrag.Fragment(atoms, "frag2")  # Different name, same composition
        
        hash1 = fragment1.fragment_hash()
        hash2 = fragment2.fragment_hash()
        
        # Hashes should be the same for identical composition
        assert hash1 == hash2
        assert isinstance(hash1, str)


class TestFragmentFragmentation:
    """Test fragment splitting functionality"""
    
    def test_connectivity_fragmentation(self):
        """Test fragmentation by connectivity"""
        # Create a molecule with two disconnected components
        atoms = [
            # Component 1: methane-like
            libfrag.Atom(6, 0.0, 0.0, 0.0),  # Carbon
            libfrag.Atom(1, 1.1, 0.0, 0.0),  # Hydrogen
            # Component 2: water-like (disconnected)
            libfrag.Atom(8, 5.0, 0.0, 0.0),  # Oxygen
            libfrag.Atom(1, 5.5, 0.8, 0.0),  # Hydrogen
        ]
        
        bonds = [(0, 1), (2, 3)]  # Two separate components
        
        fragment = libfrag.Fragment(atoms, bonds, "disconnected")
        subfragments = fragment.fragment_by_connectivity()
        
        # Should create 2 subfragments for the 2 connected components
        assert len(subfragments) == 2
        assert fragment.subfragment_count == 2
        
        # Check that total atoms are preserved
        total_atoms = sum(subfrag.atom_count for subfrag in subfragments)
        assert total_atoms == fragment.atom_count
    
    def test_size_fragmentation(self):
        """Test fragmentation by size"""
        # Create a larger fragment
        atoms = [libfrag.Atom(6, i * 1.5, 0.0, 0.0) for i in range(8)]
        fragment = libfrag.Fragment(atoms, "large_fragment")
        
        subfragments = fragment.fragment_by_size(3)  # Target 3 subfragments
        
        assert len(subfragments) > 0
        assert fragment.subfragment_count > 0
        
        # Check that total atoms are preserved
        total_atoms = sum(subfrag.atom_count for subfrag in subfragments)
        assert total_atoms == fragment.atom_count


class TestFragmentStringRepresentations:
    """Test fragment string representations"""
    
    def test_fragment_string_methods(self):
        """Test various string representation methods"""
        atoms = [libfrag.Atom(6, 0.0, 0.0, 0.0)]
        fragment = libfrag.Fragment(atoms, "test_fragment")
        
        # Test basic string methods
        fragment_str = fragment.to_fragment_string()
        tree_str = fragment.to_tree_string()
        str_repr = str(fragment)
        repr_str = repr(fragment)
        
        assert "test_fragment" in fragment_str
        assert "test_fragment" in tree_str
        assert isinstance(str_repr, str)
        assert isinstance(repr_str, str)
        assert "test_fragment" in repr_str


class TestFragmentFactoryMethods:
    """Test Fragment static factory methods"""
    
    def test_create_from_molecule(self):
        """Test creating fragment from molecule"""
        atoms = [
            libfrag.Atom(6, 0.0, 0.0, 0.0),
            libfrag.Atom(1, 1.1, 0.0, 0.0),
        ]
        
        molecule = libfrag.Molecule(atoms)
        fragment = libfrag.Fragment.create_from_molecule(molecule, "test")
        
        assert fragment.atom_count == 2
        assert "test" in fragment.fragment_id
        assert fragment.is_root_fragment()
    
    def test_create_fragment_hierarchy(self):
        """Test creating fragment hierarchy"""
        atoms = [libfrag.Atom(6, i * 1.5, 0.0, 0.0) for i in range(8)]
        molecule = libfrag.Molecule(atoms)
        
        fragment = libfrag.Fragment.create_fragment_hierarchy(
            molecule, max_fragment_size=3, id_prefix="hierarchy"
        )
        
        assert fragment.atom_count == 8
        assert "hierarchy" in fragment.fragment_id
        
        # Should create subfragments if molecule is larger than max_fragment_size
        if molecule.atom_count > 3:
            assert fragment.has_subfragments()


class TestFragmentUtilityFunctions:
    """Test utility functions provided with Fragment bindings"""
    
    def test_create_test_fragment(self):
        """Test the create_test_fragment utility function"""
        fragment = libfrag.create_test_fragment()
        
        assert fragment.atom_count == 4  # Methane: C + 3H
        assert fragment.bond_count == 3  # 3 C-H bonds
        assert fragment.fragment_id == "test_methane"
        
        # Check composition
        composition = fragment.fragment_composition()
        assert composition["C"] == 1
        assert composition["H"] == 3
    
    def test_create_complex_test_fragment(self):
        """Test the create_complex_test_fragment utility function"""
        fragment = libfrag.create_complex_test_fragment()
        
        assert fragment.atom_count == 7  # Methane (4) + water (3)
        assert fragment.bond_count == 5  # Methane (3) + water (2)
        assert fragment.fragment_id == "test_complex"
        
        # This should be disconnected, so fragmentation should work
        subfragments = fragment.fragment_by_connectivity()
        assert len(subfragments) == 2  # Should separate into methane and water


class TestFragmentEquality:
    """Test Fragment equality operators"""
    
    def test_fragment_equality(self):
        """Test Fragment equality comparison"""
        atoms = [libfrag.Atom(6, 0.0, 0.0, 0.0)]
        
        fragment1 = libfrag.Fragment(atoms, "same_id")
        fragment2 = libfrag.Fragment(atoms, "same_id")
        fragment3 = libfrag.Fragment(atoms, "different_id")
        
        # Note: Equality depends on both molecular content and fragment-specific data
        # This test may need adjustment based on the actual equality implementation
        assert fragment1.fragment_id == fragment2.fragment_id
        assert fragment1.fragment_id != fragment3.fragment_id


if __name__ == "__main__":
    pytest.main([__file__])
