#include <iostream>
#include <vector>
#include <iomanip>
#include "libfrag/atom.hpp"
#include "libfrag/bond.hpp"
#include "libfrag/molecule.hpp"

using namespace libfrag;

void demo_bond_creation() {
    std::cout << "=== Bond Creation and Properties ===" << std::endl;
    
    // Create atoms
    Atom carbon(6, 0.0, 0.0, 0.0);
    Atom hydrogen("H", 1.089, 0.0, 0.0);
    Atom nitrogen("N", 0.0, 1.5, 0.0);
    
    // Create different types of bonds
    Bond ch_bond(carbon, hydrogen, Bond::BondType::SINGLE);
    Bond cn_bond(carbon, nitrogen, 1.5);  // Partial double bond character
    
    std::cout << "C-H bond: " << ch_bond.to_string() << std::endl;
    std::cout << "C-N bond: " << cn_bond.to_string() << std::endl;
    
    // Test bond properties
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nC-H bond length: " << ch_bond.bond_length() << " Angstroms" << std::endl;
    std::cout << "C-H bond order: " << ch_bond.bond_order() << std::endl;
    std::cout << "C-H is polar: " << (ch_bond.is_polar_bond() ? "Yes" : "No") << std::endl;
    std::cout << "C-H is strong bond: " << (ch_bond.is_strong_bond() ? "Yes" : "No") << std::endl;
    
    // Set quantum properties
    ch_bond.set_property("wiberg_bond_order", 0.98);
    ch_bond.set_property("mayer_bond_order", 0.95);
    ch_bond.set_bond_strength(413.0);  // C-H bond strength in kcal/mol
    
    std::cout << "\nC-H Wiberg bond order: ";
    auto wiberg = ch_bond.get_property("wiberg_bond_order");
    if (wiberg) std::cout << *wiberg << std::endl;
    
    if (ch_bond.bond_strength()) {
        std::cout << "C-H bond strength: " << *ch_bond.bond_strength() << " kcal/mol" << std::endl;
    }
    
    std::cout << "C-H properties count: " << ch_bond.properties().size() << std::endl;
    std::cout << std::endl;
}

void demo_molecule_creation() {
    std::cout << "=== Molecule Creation ===" << std::endl;
    
    // Create empty molecule and add atoms
    Molecule methane_molecule;
    std::cout << "Empty molecule: " << methane_molecule.to_string() << std::endl;
    
    // Add atoms individually
    auto carbon_index = methane_molecule.add_atom(6, 0.0, 0.0, 0.0);
    auto h1_index = methane_molecule.add_atom("H", 1.089, 0.0, 0.0);
    auto h2_index = methane_molecule.add_atom("H", -0.363, 1.027, 0.0);
    auto h3_index = methane_molecule.add_atom("H", -0.363, -0.513, 0.889);
    auto h4_index = methane_molecule.add_atom("H", -0.363, -0.513, -0.889);
    
    std::cout << "Molecule with atoms: " << methane_molecule.to_string() << std::endl;
    std::cout << "Atom count: " << methane_molecule.atom_count() << std::endl;
    
    // Add bonds
    methane_molecule.add_bond(carbon_index, h1_index, Bond::BondType::SINGLE);
    methane_molecule.add_bond(carbon_index, h2_index, Bond::BondType::SINGLE);
    methane_molecule.add_bond(carbon_index, h3_index, Bond::BondType::SINGLE);
    methane_molecule.add_bond(carbon_index, h4_index, Bond::BondType::SINGLE);
    
    std::cout << "Methane molecule: " << methane_molecule.to_string() << std::endl;
    std::cout << "Bond count: " << methane_molecule.bond_count() << std::endl;
    std::cout << std::endl;
}

void demo_molecule_from_atoms() {
    std::cout << "=== Molecule from Atom List ===" << std::endl;
    
    // Create water molecule
    std::vector<Atom> water_atoms = {
        Atom("O", 0.0, 0.0, 0.0),
        Atom("H", 0.757, 0.586, 0.0),
        Atom("H", -0.757, 0.586, 0.0)
    };
    
    // Define connectivity (atom index pairs)
    std::vector<std::pair<std::size_t, std::size_t>> water_bonds = {{0, 1}, {0, 2}};
    std::vector<Bond::BondType> bond_types = {Bond::BondType::SINGLE, Bond::BondType::SINGLE};
    
    Molecule water_molecule(water_atoms, water_bonds, bond_types);
    std::cout << "Water molecule: " << water_molecule.to_string() << std::endl;
    std::cout << "Water is connected: " << (water_molecule.is_connected() ? "Yes" : "No") << std::endl;
    
    // Test molecular properties
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nWater molecular mass: " << water_molecule.molecular_mass() << " amu" << std::endl;
    std::cout << "Water molecular charge: " << water_molecule.molecular_charge() << std::endl;
    
    auto com = water_molecule.center_of_mass();
    std::cout << "Water center of mass: (" << com[0] << ", " << com[1] << ", " << com[2] << ")" << std::endl;
    
    // Test coordination numbers
    std::cout << "O coordination number: " << water_molecule.coordination_number(0) << std::endl;
    std::cout << "H coordination number: " << water_molecule.coordination_number(1) << std::endl;
    
    // Test bonded atoms
    auto bonded_to_oxygen = water_molecule.bonded_atoms(0);
    std::cout << "Atoms bonded to oxygen: ";
    for (auto atom_idx : bonded_to_oxygen) {
        std::cout << atom_idx << " ";
    }
    std::cout << std::endl << std::endl;
}

void demo_molecular_geometry() {
    std::cout << "=== Molecular Geometry Operations ===" << std::endl;
    
    // Create linear molecule using static method
    std::vector<std::string> element_symbols = {"C", "N", "O"};
    Molecule linear_molecule = Molecule::create_linear_molecule(element_symbols, 1.2);
    
    std::cout << "Linear CNO molecule: " << linear_molecule.to_string() << std::endl;
    
    // Test geometry properties
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Molecular diameter: " << linear_molecule.molecular_diameter() << " Angstroms" << std::endl;
    
    auto bbox = linear_molecule.bounding_box();
    std::cout << "Bounding box: [" << bbox[0] << ", " << bbox[1] << ", " << bbox[2] 
              << "] to [" << bbox[3] << ", " << bbox[4] << ", " << bbox[5] << "]" << std::endl;
    
    // Center at origin
    auto original_com = linear_molecule.center_of_mass();
    std::cout << "Original center of mass: (" << original_com[0] << ", " 
              << original_com[1] << ", " << original_com[2] << ")" << std::endl;
    
    linear_molecule.center_at_origin();
    auto new_com = linear_molecule.center_of_mass();
    std::cout << "New center of mass: (" << new_com[0] << ", " 
              << new_com[1] << ", " << new_com[2] << ")" << std::endl;
    
    // Translate molecule
    linear_molecule.translate({1.0, 2.0, 3.0});
    auto translated_com = linear_molecule.center_of_mass();
    std::cout << "Translated center of mass: (" << translated_com[0] << ", " 
              << translated_com[1] << ", " << translated_com[2] << ")" << std::endl;
    std::cout << std::endl;
}

void demo_bond_analysis() {
    std::cout << "=== Bond Analysis ===" << std::endl;
    
    // Create ethane molecule manually
    std::vector<Atom> ethane_atoms = {
        Atom("C", 0.0, 0.0, 0.0),      // 0
        Atom("C", 1.54, 0.0, 0.0),     // 1
        Atom("H", -0.51, 1.02, 0.0),   // 2
        Atom("H", -0.51, -0.51, 0.88), // 3
        Atom("H", -0.51, -0.51, -0.88),// 4
        Atom("H", 2.05, 1.02, 0.0),    // 5
        Atom("H", 2.05, -0.51, 0.88),  // 6
        Atom("H", 2.05, -0.51, -0.88)  // 7
    };
    
    Molecule ethane_molecule(ethane_atoms);
    
    // Add C-C bond
    auto cc_bond_index = ethane_molecule.add_bond(0, 1, Bond::BondType::SINGLE);
    
    // Add C-H bonds
    for (std::size_t h_index = 2; h_index < 5; ++h_index) {
        ethane_molecule.add_bond(0, h_index, Bond::BondType::SINGLE);
    }
    for (std::size_t h_index = 5; h_index < 8; ++h_index) {
        ethane_molecule.add_bond(1, h_index, Bond::BondType::SINGLE);
    }
    
    std::cout << "Ethane: " << ethane_molecule.to_string() << std::endl;
    
    // Analyze C-C bond
    Bond& cc_bond = ethane_molecule.bond(cc_bond_index);
    std::cout << "\nC-C bond: " << cc_bond.to_string() << std::endl;
    std::cout << "C-C bond length: " << std::fixed << std::setprecision(3) 
              << cc_bond.bond_length() << " Angstroms" << std::endl;
    std::cout << "C-C involves carbon: " << (cc_bond.involves_atom_type(6) ? "Yes" : "No") << std::endl;
    std::cout << "C-C involves hydrogen: " << (cc_bond.involves_atom_type(1) ? "Yes" : "No") << std::endl;
    
    // Set bond properties
    cc_bond.set_property("rotation_barrier", 2.9);  // kcal/mol
    cc_bond.set_equilibrium_bond_length(1.54);
    cc_bond.set_bond_strength(83.0);  // kcal/mol
    
    auto rotation_barrier = cc_bond.get_property("rotation_barrier");
    if (rotation_barrier) {
        std::cout << "C-C rotation barrier: " << *rotation_barrier << " kcal/mol" << std::endl;
    }
    
    std::cout << "\nDetailed C-C bond information:" << std::endl;
    std::cout << cc_bond.to_detailed_string() << std::endl;
}

void demo_qc_compatibility() {
    std::cout << "=== QC Software Compatibility ===" << std::endl;
    
    // Create a molecule for QC calculations
    Molecule qc_molecule = Molecule::create_linear_molecule({"C", "N", "O"}, 1.2);
    std::cout << "Linear CNO molecule: " << qc_molecule.to_string() << std::endl;
    
    // Get data in QC format
    auto atomic_numbers = qc_molecule.atomic_numbers();
    auto symbols = qc_molecule.element_symbols();
    
    std::cout << "Atomic numbers: ";
    for (auto z : atomic_numbers) std::cout << z << " ";
    std::cout << std::endl;
    
    std::cout << "Element symbols: ";
    for (const auto& symbol : symbols) std::cout << symbol << " ";
    std::cout << std::endl;
    
    // Set molecular properties for QC calculations
    qc_molecule.set_molecular_multiplicity(2);  // Doublet
    qc_molecule.set_property("total_energy", -168.123456);  // Hartrees
    qc_molecule.set_property("homo_energy", -0.345);
    qc_molecule.set_property("lumo_energy", 0.123);
    
    std::cout << "\nMolecular multiplicity: " << qc_molecule.molecular_multiplicity() << std::endl;
    
    auto total_energy = qc_molecule.get_property("total_energy");
    if (total_energy) {
        std::cout << "Total energy: " << *total_energy << " hartrees" << std::endl;
    }
    
    auto homo_energy = qc_molecule.get_property("homo_energy");
    if (homo_energy) {
        std::cout << "HOMO energy: " << *homo_energy << " hartrees" << std::endl;
    }
    
    // Generate XYZ format
    std::string xyz_string = qc_molecule.to_xyz_string("CNO linear molecule for QC calculation");
    std::cout << "\nXYZ format:\n" << xyz_string << std::endl;
}

void demo_molecular_fragments() {
    std::cout << "=== Molecular Fragmentation ===" << std::endl;
    
    // Create two separate fragments
    std::vector<Atom> all_atoms = {
        Atom("C", 0.0, 0.0, 0.0),   // Fragment 1
        Atom("H", 1.0, 0.0, 0.0),
        Atom("N", 5.0, 0.0, 0.0),   // Fragment 2 (far away)
        Atom("H", 6.0, 0.0, 0.0)
    };
    
    // Only bonds within fragments (no bonds between fragments)
    std::vector<std::pair<std::size_t, std::size_t>> fragment_bonds = {{0, 1}, {2, 3}};
    
    Molecule disconnected_molecule(all_atoms, fragment_bonds);
    std::cout << "Disconnected molecule: " << disconnected_molecule.to_string() << std::endl;
    std::cout << "Is connected: " << (disconnected_molecule.is_connected() ? "Yes" : "No") << std::endl;
    
    // Get separate fragments
    auto fragments = disconnected_molecule.get_fragments();
    std::cout << "Number of fragments: " << fragments.size() << std::endl;
    
    for (std::size_t i = 0; i < fragments.size(); ++i) {
        std::cout << "Fragment " << i << ": " << fragments[i].to_string() << std::endl;
        std::cout << "  Atoms: " << fragments[i].atom_count() << std::endl;
        std::cout << "  Bonds: " << fragments[i].bond_count() << std::endl;
    }
    std::cout << std::endl;
}

void demo_xyz_file_io() {
    std::cout << "=== XYZ File I/O ===" << std::endl;
    
    // Create sample XYZ string
    std::string xyz_content = R"(3
Water molecule
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000)";
    
    // Parse from XYZ
    Molecule water_from_xyz = Molecule::from_xyz_string(xyz_content);
    std::cout << "Water from XYZ: " << water_from_xyz.to_string() << std::endl;
    
    // Convert back to XYZ
    std::string xyz_output = water_from_xyz.to_xyz_string("Reconstructed water molecule");
    std::cout << "XYZ output:\n" << xyz_output << std::endl;
}

int main() {
    std::cout << "LibFrag Bond and Molecule Classes Demo" << std::endl;
    std::cout << "======================================" << std::endl << std::endl;
    
    try {
        demo_bond_creation();
        demo_molecule_creation();
        demo_molecule_from_atoms();
        demo_molecular_geometry();
        demo_bond_analysis();
        demo_qc_compatibility();
        demo_molecular_fragments();
        demo_xyz_file_io();
        
        std::cout << "======================================" << std::endl;
        std::cout << "All demos completed successfully!" << std::endl;
        std::cout << "\nThe Bond and Molecule classes are ready for use with:" << std::endl;
        std::cout << "- QCElemental for molecule data structures" << std::endl;
        std::cout << "- QCEngine for quantum chemistry calculations" << std::endl;
        std::cout << "- MolSSI Driver Interface for workflow management" << std::endl;
        std::cout << "- Fragment-based quantum chemistry methods" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
