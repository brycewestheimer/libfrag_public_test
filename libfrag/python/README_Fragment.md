# Fragment Class Python API

The Fragment class extends the Molecule class with fragment-specific functionality using the composite design pattern. This document provides a comprehensive guide to using the Fragment class from Python.

## Installation and Import

```python
import libfrag
```

## Core Classes

### Fragment

The main Fragment class that extends Molecule with hierarchical fragmentation capabilities.

```python
# Create from atoms and bonds
atoms = [libfrag.Atom(6, 0.0, 0.0, 0.0), libfrag.Atom(1, 1.1, 0.0, 0.0)]
bonds = [(0, 1)]
fragment = libfrag.Fragment(atoms, bonds, "my_fragment")
```

### FragmentLink

Represents covalent bonds between different fragments.

```python
link = libfrag.FragmentLink(
    source_fragment_id="frag1", 
    target_fragment_id="frag2",
    source_atom_index=0, 
    target_atom_index=1,
    bond_type=libfrag.BondType.SINGLE,
    bond_order=1.0
)
```

## Basic Usage

### Creating Fragments

```python
# Empty fragment
fragment = libfrag.Fragment()

# Fragment with ID
fragment = libfrag.Fragment("my_fragment")

# From atoms
atoms = [libfrag.Atom(6, 0, 0, 0), libfrag.Atom(1, 1, 0, 0)]
fragment = libfrag.Fragment(atoms, "simple_fragment")

# From atoms and bonds
bonds = [(0, 1)]
fragment = libfrag.Fragment(atoms, bonds, "bonded_fragment")

# From existing molecule
molecule = libfrag.Molecule(atoms)
fragment = libfrag.Fragment(molecule, "from_molecule")
```

### Fragment Properties

```python
# Basic properties (inherited from Molecule)
print(f"Atoms: {fragment.atom_count}")
print(f"Bonds: {fragment.bond_count}")
print(f"Is empty: {fragment.is_empty()}")

# Fragment-specific properties
print(f"Fragment ID: {fragment.fragment_id}")
print(f"Generation level: {fragment.generation_level}")
print(f"Is root: {fragment.is_root_fragment()}")
print(f"Subfragment count: {fragment.subfragment_count}")

# Calculate properties
properties = fragment.calculate_fragment_properties()
for name, value in properties.items():
    print(f"{name}: {value}")

# Get composition
composition = fragment.fragment_composition()
print(f"Composition: {composition}")  # e.g., {'C': 1, 'H': 4}
```

## Hierarchical Structure (Composite Pattern)

### Managing Subfragments

```python
parent = libfrag.Fragment("parent")
child1 = libfrag.Fragment("child1")
child2 = libfrag.Fragment("child2")

# Add subfragments
parent.add_subfragment(child1)
parent.add_subfragment(child2)

# Access subfragments
print(f"Subfragment count: {parent.subfragment_count}")
print(f"Has subfragments: {parent.has_subfragments()}")

# Get specific subfragment
first_child = parent.subfragment(0)
print(f"First child ID: {first_child.fragment_id}")

# Get all subfragments
all_children = parent.subfragments
for child in all_children:
    print(f"Child: {child.fragment_id}")

# Remove subfragments
parent.remove_subfragment(0)
parent.clear_subfragments()
```

### Parent-Child Relationships

```python
# Check parent relationship
if not child1.is_root_fragment():
    parent = child1.parent_fragment
    if parent:
        print(f"Parent: {parent.fragment_id}")

# Get root fragment
root = child1.get_root_fragment()
print(f"Root: {root.fragment_id}")

# Total atoms in entire tree
total_atoms = parent.total_atom_count()
```

## Fragment Splitting

### By Connectivity

Split disconnected molecular components into separate fragments:

```python
# Create molecule with disconnected parts
atoms = [
    # Methane
    libfrag.Atom(6, 0.0, 0.0, 0.0),
    libfrag.Atom(1, 1.1, 0.0, 0.0),
    # Water (disconnected)
    libfrag.Atom(8, 5.0, 0.0, 0.0),
    libfrag.Atom(1, 5.5, 0.8, 0.0),
]
bonds = [(0, 1), (2, 3)]

fragment = libfrag.Fragment(atoms, bonds, "mixed")
subfragments = fragment.fragment_by_connectivity()

print(f"Split into {len(subfragments)} components")
for i, subfrag in enumerate(subfragments):
    print(f"  Component {i}: {subfrag.atom_count} atoms")
```

### By Size

Split fragment into roughly equal-sized pieces:

```python
# Create large fragment
atoms = [libfrag.Atom(6, i * 1.5, 0, 0) for i in range(10)]
fragment = libfrag.Fragment(atoms, "large_fragment")

# Split into 3 pieces
subfragments = fragment.fragment_by_size(3)
print(f"Split into {len(subfragments)} size-based fragments")
```

### By Bond Breaking

Split by breaking specific bonds:

```python
# Fragment with bonds
atoms = [libfrag.Atom(6, i, 0, 0) for i in range(4)]
bonds = [(0, 1), (1, 2), (2, 3)]  # Linear chain
fragment = libfrag.Fragment(atoms, bonds, "chain")

# Break middle bond
subfragments = fragment.fragment_by_bond_breaking([1])  # Break bond index 1
print(f"Bond breaking created {len(subfragments)} fragments")
```

### Custom Fragmentation

Use a custom function to define fragmentation logic:

```python
def custom_fragmenter(fragment):
    """Group atoms by element type"""
    carbon_atoms = []
    hydrogen_atoms = []
    
    for i in range(fragment.atom_count):
        atom = fragment.atom(i)
        if atom.atomic_number == 6:
            carbon_atoms.append(i)
        elif atom.atomic_number == 1:
            hydrogen_atoms.append(i)
    
    groups = []
    if carbon_atoms:
        groups.append(carbon_atoms)
    if hydrogen_atoms:
        groups.append(hydrogen_atoms)
    
    return groups

subfragments = fragment.fragment_by_custom_function(custom_fragmenter)
```

## Inter-Fragment Links

### Managing Links Between Fragments

```python
fragment1 = libfrag.Fragment("fragment1")
fragment2 = libfrag.Fragment("fragment2")

# Create link
link = libfrag.FragmentLink(
    "fragment1", "fragment2",
    0, 0,  # Connect atom 0 in each fragment
    libfrag.BondType.SINGLE,
    1.0
)

# Add link to fragment
fragment1.add_fragment_link(link)

# Query links
print(f"Is linked to fragment2: {fragment1.is_linked_to('fragment2')}")
print(f"Linked fragments: {list(fragment1.linked_fragment_ids)}")

# Find specific links
link_indices = fragment1.find_links_to_fragment("fragment2")
print(f"Link indices to fragment2: {link_indices}")

# Access all links
for i, link in enumerate(fragment1.fragment_links):
    print(f"Link {i}: {link}")
```

## Factory Methods

### Create from Molecule

```python
molecule = libfrag.Molecule(atoms)

# Simple creation
fragment = libfrag.Fragment.create_from_molecule(molecule, "test")

# Hierarchical creation with automatic fragmentation
hierarchy = libfrag.Fragment.create_fragment_hierarchy(
    molecule, 
    max_fragment_size=5,  # Max atoms per fragment
    id_prefix="auto"
)
```

## Utility Functions

### Pre-built Test Fragments

```python
# Simple methane fragment for testing
methane = libfrag.create_test_fragment()
print(f"Test fragment: {methane.fragment_id}, atoms: {methane.atom_count}")

# Complex fragment with multiple components
complex_frag = libfrag.create_complex_test_fragment()
print(f"Complex fragment: {complex_frag.fragment_id}, atoms: {complex_frag.atom_count}")

# Fragment the complex molecule
subfragments = complex_frag.fragment_by_connectivity()
print(f"Split into {len(subfragments)} components")
```

## String Representations

### Display Fragment Information

```python
# Basic string representation
print(str(fragment))

# Detailed fragment information
print(fragment.to_fragment_string())

# Tree structure representation
print(fragment.to_tree_string())

# Fragment hash for comparison
hash_str = fragment.fragment_hash()
print(f"Fragment hash: {hash_str}")
```

## Fragment Analysis

### Properties and Composition

```python
# Calculate all fragment properties
properties = fragment.calculate_fragment_properties()
important_props = ['molecular_mass', 'atom_count', 'bond_count', 'subfragment_count']
for prop in important_props:
    if prop in properties:
        print(f"{prop}: {properties[prop]}")

# Element composition
composition = fragment.fragment_composition()
molecular_formula = "".join([f"{elem}{count}" for elem, count in sorted(composition.items())])
print(f"Molecular formula: {molecular_formula}")

# Molecular properties (inherited from Molecule)
print(f"Molecular mass: {fragment.molecular_mass:.2f} amu")
print(f"Center of mass: {fragment.center_of_mass()}")
print(f"Is connected: {fragment.is_connected()}")
```

## Complete Example

```python
import libfrag

# Create a complex molecule (methane + water, disconnected)
atoms = [
    # Methane
    libfrag.Atom(6, 0.0, 0.0, 0.0),  # Carbon
    libfrag.Atom(1, 1.1, 0.0, 0.0),  # H1
    libfrag.Atom(1, -0.5, 0.9, 0.0), # H2
    libfrag.Atom(1, -0.5, -0.9, 0.0),# H3
    # Water
    libfrag.Atom(8, 5.0, 0.0, 0.0),  # Oxygen
    libfrag.Atom(1, 5.5, 0.8, 0.0),  # H4
    libfrag.Atom(1, 5.5, -0.8, 0.0), # H5
]

bonds = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6)]

# Create root fragment
root = libfrag.Fragment(atoms, bonds, "methane_water")
print(f"Created: {root.fragment_id} with {root.atom_count} atoms")

# Fragment by connectivity
components = root.fragment_by_connectivity()
print(f"Split into {len(components)} components:")

for i, comp in enumerate(components):
    composition = comp.fragment_composition()
    formula = "".join([f"{elem}{count}" for elem, count in sorted(composition.items())])
    print(f"  Component {i}: {formula} ({comp.atom_count} atoms)")

# Display tree structure
print("\\nFragment tree:")
print(root.to_tree_string())

# Create links between components
if len(components) >= 2:
    link = libfrag.FragmentLink(
        components[0].fragment_id, components[1].fragment_id,
        0, 0, libfrag.BondType.SINGLE
    )
    components[0].add_fragment_link(link)
    print(f"\\nLinked {components[0].fragment_id} to {components[1].fragment_id}")
```

## Error Handling

```python
try:
    # Fragment operations that might fail
    subfragment = fragment.subfragment(999)  # Index out of range
except IndexError as e:
    print(f"Index error: {e}")

try:
    # Invalid bond type
    link = libfrag.FragmentLink("frag1", "frag2", 0, 1, "invalid", 1.0)
except (ValueError, TypeError) as e:
    print(f"Invalid link: {e}")
```

This API provides a powerful and flexible way to work with molecular fragments in Python, supporting hierarchical organization, multiple fragmentation strategies, and inter-fragment connectivity management.
