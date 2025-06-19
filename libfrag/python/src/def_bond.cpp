#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

#include "libfrag/bond.hpp"
#include "libfrag/atom.hpp"

namespace py = pybind11;

namespace libfrag {

    void def_bond(py::module& m) {
        // Bond type enumeration
        py::enum_<Bond::BondType>(m, "BondType", R"pbdoc(
            Enumeration of chemical bond types.
            
            Used to classify different kinds of chemical bonds based on their
            electronic structure and bonding mechanism.
        )pbdoc")
            .value("SINGLE", Bond::BondType::SINGLE, "Single covalent bond (bond order = 1)")
            .value("DOUBLE", Bond::BondType::DOUBLE, "Double covalent bond (bond order = 2)")
            .value("TRIPLE", Bond::BondType::TRIPLE, "Triple covalent bond (bond order = 3)")
            .value("AROMATIC", Bond::BondType::AROMATIC, "Aromatic bond (delocalized electrons)")
            .value("COORDINATE", Bond::BondType::COORDINATE, "Coordinate/dative bond")
            .value("HYDROGEN", Bond::BondType::HYDROGEN, "Hydrogen bond (weak interaction)")
            .value("VAN_DER_WAALS", Bond::BondType::VAN_DER_WAALS, "van der Waals interaction")
            .value("IONIC", Bond::BondType::IONIC, "Ionic bond")
            .value("METALLIC", Bond::BondType::METALLIC, "Metallic bond")
            .value("CUSTOM", Bond::BondType::CUSTOM, "User-defined bond type")
            .export_values();

        // Bond class
        py::class_<Bond>(m, "Bond", R"pbdoc(
            Chemical bond class representing connection between two atoms.
            
            The Bond class represents a chemical bond between two atoms in a molecular system.
            It stores references to the bonded atoms, bond properties such as order and length,
            and additional quantum mechanical properties. Designed for quantum chemistry
            applications and compatibility with MolSSI tools.
            
            Examples
            --------
            Create bonds between atoms:
            
            >>> atom1 = Atom(6, 0.0, 0.0, 0.0)  # Carbon
            >>> atom2 = Atom(1, 1.089, 0.0, 0.0)  # Hydrogen
            >>> bond = Bond(atom1, atom2, BondType.SINGLE)
            
            >>> # Check bond properties
            >>> print(bond.bond_length())  # Distance between atoms
            >>> print(bond.bond_order)     # 1.0 for single bond
            >>> print(bond.is_polar_bond())  # True for C-H
            
            >>> # Set quantum properties
            >>> bond.set_property("wiberg_bond_order", 0.98)
            >>> bond.set_property("mayer_bond_order", 0.95)
        )pbdoc")
        
        // Constructors
        .def(py::init<const Atom&, const Atom&, Bond::BondType>(),
             py::arg("first_atom"), py::arg("second_atom"), 
             py::arg("bond_type") = Bond::BondType::SINGLE,
             "Create bond between two atoms with specified type")
        .def(py::init<const Atom&, const Atom&, double>(),
             py::arg("first_atom"), py::arg("second_atom"), py::arg("bond_order"),
             "Create bond with explicit bond order")
        
        // Atom access methods
        .def_property_readonly("first_atom", &Bond::first_atom,
                              "First atom in the bond", py::return_value_policy::reference)
        .def_property_readonly("second_atom", &Bond::second_atom,
                              "Second atom in the bond", py::return_value_policy::reference)
        .def("get_other_atom", &Bond::get_other_atom, py::arg("reference_atom"),
             "Get the other atom in bond given one atom",
             py::return_value_policy::reference)
        .def("contains_atom", &Bond::contains_atom, py::arg("atom"),
             "Check if bond contains specified atom")
        
        // Bond properties
        .def_property("bond_type", &Bond::bond_type, &Bond::set_bond_type,
                     "Type of chemical bond")
        .def_property("bond_order", &Bond::bond_order, &Bond::set_bond_order,
                     "Bond order (fractional values allowed)")
        .def_property_readonly("bond_length", &Bond::bond_length,
                              "Current bond length from atomic coordinates (Angstroms)")
        .def_property("equilibrium_bond_length", &Bond::equilibrium_bond_length, 
                     &Bond::set_equilibrium_bond_length,
                     "Equilibrium bond length (Angstroms)")
        .def_property("bond_strength", &Bond::bond_strength, &Bond::set_bond_strength,
                     "Bond dissociation energy (kcal/mol)")
        
        // Quantum properties
        .def_property_readonly("properties", &Bond::properties,
                              "Dictionary of quantum mechanical properties")
        .def("get_property", &Bond::get_property, py::arg("property_name"),
             "Get quantum property value by name")
        .def("set_property", &Bond::set_property, 
             py::arg("property_name"), py::arg("property_value"),
             "Set quantum mechanical property")
        
        // Bond classification methods
        .def("is_strong_bond", &Bond::is_strong_bond,
             "Check if bond is strong (covalent/ionic)")
        .def("is_weak_bond", &Bond::is_weak_bond,
             "Check if bond is weak (hydrogen/van der Waals)")
        .def("involves_atom_type", &Bond::involves_atom_type,
             py::arg("atomic_number"),
             "Check if bond involves specified atom type")
        .def("is_polar_bond", &Bond::is_polar_bond,
             "Check if bond is polar (different electronegativity)")
        
        // String representations
        .def("to_string", &Bond::to_string,
             "Brief string representation of bond")
        .def("to_detailed_string", &Bond::to_detailed_string,
             "Detailed string with all bond properties")
        .def("__str__", &Bond::to_string)
        .def("__repr__", [](const Bond& bond) {
             return "Bond(" + bond.first_atom().symbol() + "-" + 
                    bond.second_atom().symbol() + ", " + 
                    Bond::bond_type_to_string(bond.bond_type()) + ")";
        })
        
        // Static utility methods
        .def_static("bond_type_to_string", &Bond::bond_type_to_string,
                   py::arg("bond_type"), "Convert bond type enum to string")
        .def_static("string_to_bond_type", &Bond::string_to_bond_type,
                   py::arg("type_string"), "Convert string to bond type enum")
        .def_static("get_typical_bond_order", &Bond::get_typical_bond_order,
                   py::arg("bond_type"), "Get typical bond order for bond type")
        
        // Operators
        .def(py::self == py::self)
        .def(py::self != py::self);
    }

} // namespace libfrag
