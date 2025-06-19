#include "libfrag/fragment.hpp"
#include "libfrag/atom.hpp"
#include <iostream>
#include <vector>

int main() {
    using namespace libfrag;
    
    // Create some atoms for testing
    std::vector<Atom> atoms = {
        Atom(6, 0.0, 0.0, 0.0),    // Carbon
        Atom(1, 1.1, 0.0, 0.0),    // Hydrogen
        Atom(1, -0.5, 0.9, 0.0),   // Hydrogen
        Atom(1, -0.5, -0.9, 0.0),  // Hydrogen
        Atom(8, 2.5, 0.0, 0.0),    // Oxygen (separate fragment)
        Atom(1, 3.0, 0.8, 0.0)     // Hydrogen
    };
    
    // Create bond connectivity for methane + water
    std::vector<std::pair<std::size_t, std::size_t>> bonds = {
        {0, 1}, {0, 2}, {0, 3},  // Methane bonds
        {4, 5}                   // Water bond
    };
    
    std::cout << "=== Fragment Class Demo ===" << std::endl;
    
    // Create root fragment
    auto root_fragment = std::make_shared<Fragment>(atoms, bonds, "root_molecule");
    std::cout << "Created root fragment with " << root_fragment->atom_count() 
              << " atoms and " << root_fragment->bond_count() << " bonds" << std::endl;
    
    // Test fragment properties
    auto properties = root_fragment->calculate_fragment_properties();
    std::cout << "\nFragment Properties:" << std::endl;
    for (const auto& [name, value] : properties) {
        std::cout << "  " << name << ": " << value << std::endl;
    }
    
    // Test composition
    auto composition = root_fragment->fragment_composition();
    std::cout << "\nFragment Composition:" << std::endl;
    for (const auto& [element, count] : composition) {
        std::cout << "  " << element << ": " << count << std::endl;
    }
    
    // Test fragmentation by connectivity
    std::cout << "\n=== Fragmentation by Connectivity ===" << std::endl;
    auto subfragments = root_fragment->fragment_by_connectivity();
    std::cout << "Created " << subfragments.size() << " subfragments:" << std::endl;
    
    for (std::size_t i = 0; i < subfragments.size(); ++i) {
        const auto& subfrag = subfragments[i];
        std::cout << "  Subfragment " << i << " (ID: " << subfrag->fragment_id() 
                  << "): " << subfrag->atom_count() << " atoms, " 
                  << subfrag->bond_count() << " bonds" << std::endl;
        
        auto sub_composition = subfrag->fragment_composition();
        std::cout << "    Composition: ";
        for (const auto& [element, count] : sub_composition) {
            std::cout << element << count << " ";
        }
        std::cout << std::endl;
    }
    
    // Test fragment tree representation
    std::cout << "\n=== Fragment Tree ===" << std::endl;
    std::cout << root_fragment->to_tree_string() << std::endl;
    
    // Test fragment links
    std::cout << "=== Fragment Links Demo ===" << std::endl;
    if (subfragments.size() >= 2) {
        // Create a link between first two subfragments
        FragmentLink link(subfragments[0]->fragment_id(), 
                         subfragments[1]->fragment_id(),
                         0, 0,  // Connect first atoms
                         Bond::BondType::SINGLE, 1.5);
        
        subfragments[0]->add_fragment_link(link);
        std::cout << "Added link: " << link.to_string() << std::endl;
        std::cout << "Fragment " << subfragments[0]->fragment_id() 
                  << " is linked to: ";
        auto linked_ids = subfragments[0]->linked_fragment_ids();
        for (const auto& id : linked_ids) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
    }
    
    // Test static factory methods
    std::cout << "\n=== Static Factory Methods ===" << std::endl;
    
    // Create a simple molecule for testing
    Molecule simple_molecule(atoms);
    auto fragment_from_molecule = Fragment::create_from_molecule(simple_molecule, "test");
    std::cout << "Created fragment from molecule: " << fragment_from_molecule->fragment_id() 
              << " with " << fragment_from_molecule->atom_count() << " atoms" << std::endl;
    
    // Test hierarchical fragmentation
    auto hierarchical_fragment = Fragment::create_fragment_hierarchy(simple_molecule, 3, "hierarchy");
    std::cout << "Created hierarchical fragment: " << hierarchical_fragment->fragment_id() 
              << " with " << hierarchical_fragment->subfragment_count() << " subfragments" << std::endl;
    
    std::cout << "\nHierarchical fragment tree:" << std::endl;
    std::cout << hierarchical_fragment->to_tree_string() << std::endl;
    
    std::cout << "\n=== Fragment Demo Complete ===" << std::endl;
    
    return 0;
}
