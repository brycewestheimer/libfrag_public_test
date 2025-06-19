#include "libfrag/fragment.hpp"
#include "libfrag/atom.hpp"
#include <doctest/doctest.h>
#include <vector>

using namespace libfrag;

TEST_CASE("Fragment Basic Construction") {
    SUBCASE("Default constructor") {
        Fragment fragment;
        CHECK(fragment.is_empty());
        CHECK(fragment.fragment_id().empty());
        CHECK(fragment.is_root_fragment());
        CHECK(fragment.generation_level() == 0);
    }
    
    SUBCASE("Constructor with ID") {
        Fragment fragment("test_fragment");
        CHECK(fragment.is_empty());
        CHECK(fragment.fragment_id() == "test_fragment");
        CHECK(fragment.is_root_fragment());
    }
    
    SUBCASE("Constructor from atoms") {
        std::vector<Atom> atoms = {
            Atom(6, 0.0, 0.0, 0.0),
            Atom(1, 1.1, 0.0, 0.0)
        };
        Fragment fragment(atoms, "methane_part");
        CHECK(fragment.atom_count() == 2);
        CHECK(fragment.fragment_id() == "methane_part");
    }
}

TEST_CASE("Fragment Molecular Properties") {
    std::vector<Atom> atoms = {
        Atom(6, 0.0, 0.0, 0.0),    // Carbon
        Atom(1, 1.1, 0.0, 0.0),    // Hydrogen
        Atom(1, -0.5, 0.9, 0.0),   // Hydrogen
        Atom(1, -0.5, -0.9, 0.0)   // Hydrogen
    };
    
    std::vector<std::pair<std::size_t, std::size_t>> bonds = {
        {0, 1}, {0, 2}, {0, 3}
    };
    
    Fragment fragment(atoms, bonds, "methane");
    
    SUBCASE("Basic properties") {
        CHECK(fragment.atom_count() == 4);
        CHECK(fragment.bond_count() == 3);
        CHECK(fragment.fragment_id() == "methane");
    }
    
    SUBCASE("Fragment composition") {
        auto composition = fragment.fragment_composition();
        CHECK(composition["C"] == 1);
        CHECK(composition["H"] == 3);
        CHECK(composition.size() == 2);
    }
    
    SUBCASE("Fragment properties calculation") {
        auto properties = fragment.calculate_fragment_properties();
        CHECK(properties["atom_count"] == 4.0);
        CHECK(properties["bond_count"] == 3.0);
        CHECK(properties["subfragment_count"] == 0.0);
        CHECK(properties["generation_level"] == 0.0);
    }
}

TEST_CASE("Fragment Hierarchy Management") {
    // Create parent fragment
    std::vector<Atom> atoms = {
        Atom(6, 0.0, 0.0, 0.0),
        Atom(1, 1.1, 0.0, 0.0),
        Atom(8, 3.0, 0.0, 0.0),
        Atom(1, 3.5, 0.8, 0.0)
    };
    
    auto parent = std::make_shared<Fragment>(atoms, "parent");
    
    SUBCASE("Subfragment management") {
        CHECK(parent->subfragment_count() == 0);
        CHECK_FALSE(parent->has_subfragments());
        
        // Create and add subfragments
        auto child1 = std::make_shared<Fragment>("child1");
        auto child2 = std::make_shared<Fragment>("child2");
        
        parent->add_subfragment(child1);
        parent->add_subfragment(child2);
        
        CHECK(parent->subfragment_count() == 2);
        CHECK(parent->has_subfragments());
        CHECK(parent->subfragment(0)->fragment_id() == "child1");
        CHECK(parent->subfragment(1)->fragment_id() == "child2");
    }
    
    SUBCASE("Parent-child relationships") {
        auto child = std::make_shared<Fragment>("child");
        parent->add_subfragment(child);
        
        CHECK(child->generation_level() == 1);
        CHECK_FALSE(child->is_root_fragment());
        
        auto child_parent = child->parent_fragment().lock();
        CHECK(child_parent != nullptr);
        CHECK(child_parent->fragment_id() == "parent");
    }
    
    SUBCASE("Subfragment removal") {
        auto child1 = std::make_shared<Fragment>("child1");
        auto child2 = std::make_shared<Fragment>("child2");
        
        parent->add_subfragment(child1);
        parent->add_subfragment(child2);
        
        CHECK(parent->subfragment_count() == 2);
        
        parent->remove_subfragment(0);
        CHECK(parent->subfragment_count() == 1);
        CHECK(parent->subfragment(0)->fragment_id() == "child2");
        
        parent->clear_subfragments();
        CHECK(parent->subfragment_count() == 0);
    }
}

TEST_CASE("Fragment Links") {
    SUBCASE("FragmentLink creation and properties") {
        FragmentLink link("frag1", "frag2", 0, 1, Bond::BondType::DOUBLE, 2.0);
        
        CHECK(link.source_fragment_id() == "frag1");
        CHECK(link.target_fragment_id() == "frag2");
        CHECK(link.source_atom_index() == 0);
        CHECK(link.target_atom_index() == 1);
        CHECK(link.bond_type() == Bond::BondType::DOUBLE);
        CHECK(link.bond_order() == 2.0);
    }
    
    SUBCASE("Fragment link management") {
        Fragment fragment("test_fragment");
        
        FragmentLink link1("test_fragment", "other1", 0, 0);
        FragmentLink link2("test_fragment", "other2", 1, 1);
        
        fragment.add_fragment_link(link1);
        fragment.add_fragment_link(link2);
        
        CHECK(fragment.fragment_links().size() == 2);
        CHECK(fragment.is_linked_to("other1"));
        CHECK(fragment.is_linked_to("other2"));
        CHECK_FALSE(fragment.is_linked_to("other3"));
        
        auto linked_ids = fragment.linked_fragment_ids();
        CHECK(linked_ids.size() == 2);
        CHECK(linked_ids.count("other1") > 0);
        CHECK(linked_ids.count("other2") > 0);
    }
}

TEST_CASE("Fragment Connectivity Analysis") {
    // Create a molecule with two disconnected components
    std::vector<Atom> atoms = {
        // Component 1: methane
        Atom(6, 0.0, 0.0, 0.0),    // 0: Carbon
        Atom(1, 1.1, 0.0, 0.0),    // 1: Hydrogen
        Atom(1, -0.5, 0.9, 0.0),   // 2: Hydrogen
        // Component 2: water (disconnected)
        Atom(8, 5.0, 0.0, 0.0),    // 3: Oxygen
        Atom(1, 5.5, 0.8, 0.0),    // 4: Hydrogen
    };
    
    std::vector<std::pair<std::size_t, std::size_t>> bonds = {
        {0, 1}, {0, 2},  // Methane bonds
        {3, 4}           // Water bond
    };
    
    Fragment fragment(atoms, bonds, "disconnected");
    
    SUBCASE("Fragmentation by connectivity") {
        auto subfragments = fragment.fragment_by_connectivity();
        
        // Should create 2 subfragments for the 2 connected components
        CHECK(subfragments.size() == 2);
        CHECK(fragment.subfragment_count() == 2);
        
        // Check that atoms are properly distributed
        std::size_t total_atoms = 0;
        for (const auto& subfrag : subfragments) {
            total_atoms += subfrag->atom_count();
            CHECK(subfrag->generation_level() == 1);
        }
        CHECK(total_atoms == fragment.atom_count());
    }
}

TEST_CASE("Fragment Static Factory Methods") {
    SUBCASE("Create from molecule") {
        std::vector<Atom> atoms = {
            Atom(6, 0.0, 0.0, 0.0),
            Atom(1, 1.1, 0.0, 0.0)
        };
        
        Molecule molecule(atoms);
        auto fragment = Fragment::create_from_molecule(molecule, "test");
        
        CHECK(fragment->atom_count() == 2);
        CHECK(fragment->fragment_id().find("test") != std::string::npos);
        CHECK(fragment->is_root_fragment());
    }
    
    SUBCASE("Create fragment hierarchy") {
        std::vector<Atom> atoms;
        // Create a larger molecule for hierarchy testing
        for (int i = 0; i < 8; ++i) {
            atoms.emplace_back(6, i * 1.5, 0.0, 0.0);  // Carbon chain
        }
        
        Molecule molecule(atoms);
        auto fragment = Fragment::create_fragment_hierarchy(molecule, 3, "hierarchy");
        
        CHECK(fragment->atom_count() == 8);
        CHECK(fragment->fragment_id().find("hierarchy") != std::string::npos);
        
        // Should create subfragments if molecule is larger than max_fragment_size
        if (molecule.atom_count() > 3) {
            CHECK(fragment->has_subfragments());
        }
    }
}

TEST_CASE("Fragment String Representations") {
    std::vector<Atom> atoms = {
        Atom(6, 0.0, 0.0, 0.0),
        Atom(1, 1.1, 0.0, 0.0)
    };
    
    Fragment fragment(atoms, "test_fragment");
    
    SUBCASE("Fragment string representation") {
        std::string fragment_str = fragment.to_fragment_string();
        CHECK(fragment_str.find("test_fragment") != std::string::npos);
        CHECK(fragment_str.find("Atoms: 2") != std::string::npos);
    }
    
    SUBCASE("Tree string representation") {
        std::string tree_str = fragment.to_tree_string();
        CHECK(tree_str.find("test_fragment") != std::string::npos);
        CHECK(tree_str.find("atoms: 2") != std::string::npos);
    }
    
    SUBCASE("Fragment hash") {
        std::string hash1 = fragment.fragment_hash();
        
        // Create identical fragment
        Fragment fragment2(atoms, "different_name");
        std::string hash2 = fragment2.fragment_hash();
        
        // Hashes should be the same for identical composition
        CHECK(hash1 == hash2);
    }
}
