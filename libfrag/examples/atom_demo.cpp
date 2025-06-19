#include <iostream>
#include <vector>
#include <iomanip>
#include "libfrag/atom.hpp"

using namespace libfrag;

void demo_basic_usage() {
    std::cout << "=== Basic Atom Usage ===" << std::endl;
    
    // Create atoms different ways
    Atom carbon(6, 0.0, 0.0, 0.0);
    Atom hydrogen("H", 1.089, 0.0, 0.0);
    
    std::cout << "Carbon: " << carbon.to_string() << std::endl;
    std::cout << "Hydrogen: " << hydrogen.to_string() << std::endl;
    
    std::cout << "Carbon atomic mass: " << carbon.atomic_mass() << " amu" << std::endl;
    std::cout << "C-H distance: " << carbon.distance_to(hydrogen) << " Angstrom" << std::endl;
    std::cout << std::endl;
}

void demo_water_molecule() {
    std::cout << "=== Water Molecule Demo ===" << std::endl;
    
    // Create water molecule
    Atom oxygen("O", 0.0, 0.0, 0.0);
    Atom hydrogen1("H", 0.757, 0.586, 0.0);
    Atom hydrogen2("H", -0.757, 0.586, 0.0);
    
    std::cout << "Water molecule coordinates:" << std::endl;
    std::cout << oxygen.to_xyz_string() << std::endl;
    std::cout << hydrogen1.to_xyz_string() << std::endl;
    std::cout << hydrogen2.to_xyz_string() << std::endl;
    
    // Calculate bond distances
    double oh1_distance = oxygen.distance_to(hydrogen1);
    double oh2_distance = oxygen.distance_to(hydrogen2);
    double hh_distance = hydrogen1.distance_to(hydrogen2);
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nBond distances:" << std::endl;
    std::cout << "O-H1: " << oh1_distance << " Angstrom" << std::endl;
    std::cout << "O-H2: " << oh2_distance << " Angstrom" << std::endl;
    std::cout << "H1-H2: " << hh_distance << " Angstrom" << std::endl;
    std::cout << std::endl;
}

void demo_quantum_properties() {
    std::cout << "=== Quantum Properties Demo ===" << std::endl;
    
    Atom carbon("C", 0.0, 0.0, 0.0);
    
    // Set various quantum properties
    carbon.set_property("mulliken_charge", -0.123);
    carbon.set_property("esp_charge", -0.145);
    carbon.set_property("spin_density", 0.0);
    carbon.set_formal_charge(0.0);
    carbon.set_partial_charge(-0.123);
    carbon.set_multiplicity(1);
    
    std::cout << "Carbon atom properties:" << std::endl;
    std::cout << "Formal charge: " << carbon.formal_charge() << std::endl;
    std::cout << "Partial charge: " << carbon.partial_charge() << std::endl;
    std::cout << "Multiplicity: " << carbon.multiplicity() << std::endl;
    
    auto mulliken = carbon.get_property("mulliken_charge");
    if (mulliken) {
        std::cout << "Mulliken charge: " << *mulliken << std::endl;
    }
    
    auto esp = carbon.get_property("esp_charge");
    if (esp) {
        std::cout << "ESP charge: " << *esp << std::endl;
    }
    
    std::cout << "Total properties stored: " << carbon.properties().size() << std::endl;
    std::cout << std::endl;
}

void demo_methane_molecule() {
    std::cout << "=== Methane Molecule Demo ===" << std::endl;
    
    // Create methane molecule
    std::vector<Atom> atoms = {
        Atom("C", 0.0, 0.0, 0.0),
        Atom("H", 1.089, 0.0, 0.0),
        Atom("H", -0.363, 1.027, 0.0),
        Atom("H", -0.363, -0.513, 0.889),
        Atom("H", -0.363, -0.513, -0.889)
    };
    
    std::cout << "Methane molecule (CH4):" << std::endl;
    for (size_t i = 0; i < atoms.size(); ++i) {
        std::cout << i << ": " << atoms[i].to_xyz_string() << std::endl;
    }
    
    // Calculate all C-H distances
    std::cout << "\nC-H bond distances:" << std::endl;
    for (size_t i = 1; i < atoms.size(); ++i) {
        double distance = atoms[0].distance_to(atoms[i]);
        std::cout << "C-H" << i << ": " << std::setprecision(3) << distance << " Angstrom" << std::endl;
    }
    
    // Calculate molecular center of mass
    double total_mass = 0.0;
    std::array<double, 3> com = {0.0, 0.0, 0.0};
    
    for (const auto& atom : atoms) {
        double mass = atom.atomic_mass();
        total_mass += mass;
        com[0] += mass * atom.x();
        com[1] += mass * atom.y();
        com[2] += mass * atom.z();
    }
    
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;
    
    std::cout << "\nMolecular properties:" << std::endl;
    std::cout << "Total mass: " << total_mass << " amu" << std::endl;
    std::cout << "Center of mass: (" << com[0] << ", " << com[1] << ", " << com[2] << ")" << std::endl;
    std::cout << std::endl;
}

void demo_static_methods() {
    std::cout << "=== Static Methods Demo ===" << std::endl;
    
    // Symbol/atomic number conversion
    std::cout << "'C' -> " << Atom::symbol_to_atomic_number("C") << std::endl;
    std::cout << "6 -> '" << Atom::atomic_number_to_symbol(6) << "'" << std::endl;
    
    // Validation
    std::cout << "Is 'C' valid? " << (Atom::is_valid_element_symbol("C") ? "Yes" : "No") << std::endl;
    std::cout << "Is 'Xx' valid? " << (Atom::is_valid_element_symbol("Xx") ? "Yes" : "No") << std::endl;
    std::cout << "Is Z=6 valid? " << (Atom::is_valid_atomic_number(6) ? "Yes" : "No") << std::endl;
    std::cout << "Is Z=200 valid? " << (Atom::is_valid_atomic_number(200) ? "Yes" : "No") << std::endl;
    std::cout << std::endl;
}

int main() {
    std::cout << "LibFrag Atom Class Demo" << std::endl;
    std::cout << "======================" << std::endl << std::endl;
    
    try {
        demo_basic_usage();
        demo_water_molecule();
        demo_quantum_properties();
        demo_methane_molecule();
        demo_static_methods();
        
        std::cout << "All demos completed successfully!" << std::endl;
        std::cout << "\nThe Atom class is ready for use with:" << std::endl;
        std::cout << "- MolSSI Driver Interface" << std::endl;
        std::cout << "- QCEngine for quantum chemistry calculations" << std::endl;
        std::cout << "- QCElemental for molecule data structures" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
