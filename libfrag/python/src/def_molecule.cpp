#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/operators.h"

#include "xtensor-python/pyarray.hpp"

#include "libfrag/molecule.hpp"
#include "libfrag/bond.hpp"
#include "libfrag/atom.hpp"

namespace py = pybind11;

namespace libfrag {

    void def_molecule(py::module& m) {
        py::class_<Molecule>(m, "Molecule", R"pbdoc(
            Comprehensive molecular system class for quantum chemistry.
            
            The Molecule class represents a complete molecular system containing atoms and bonds.
            It provides methods for molecular property calculation, geometry manipulation,
            and compatibility with quantum chemistry software packages. The class maintains
            unique atom instances and manages bonds between them using references.
            
            Key features:
            - Unique atom management with automatic indexing
            - Bond connectivity tracking
            - Molecular property calculation (center of mass, moments of inertia, etc.)
            - Quantum chemistry software compatibility (QCElemental, QCEngine)
            - Molecular geometry operations (translation, rotation)
            - Fragment analysis and manipulation
            
            Examples
            --------
            Create molecules in different ways:
            
            >>> # Empty molecule
            >>> mol = Molecule()
            
            >>> # From list of atoms
            >>> atoms = [Atom(6, 0.0, 0.0, 0.0), Atom(1, 1.089, 0.0, 0.0)]
            >>> mol = Molecule(atoms)
            
            >>> # Add bonds
            >>> bond_idx = mol.add_bond(0, 1, BondType.SINGLE)
            
            >>> # Molecular properties
            >>> print(f"Mass: {mol.molecular_mass()} amu")
            >>> print(f"Charge: {mol.molecular_charge()}")
            >>> print(f"Center of mass: {mol.center_of_mass()}")
            
            >>> # QC software compatibility
            >>> geometry = mol.geometry_array()  # For QCElemental
            >>> symbols = mol.element_symbols()  # For QCEngine
        )pbdoc")
        
        // Constructors
        .def(py::init<>(), "Create empty molecule")
        .def(py::init<const std::vector<Atom>&>(), py::arg("atoms"),
             "Create molecule from list of atoms (no bonds)")
        .def(py::init<const std::vector<Atom>&, 
                     const std::vector<std::pair<std::size_t, std::size_t>>&,
                     const std::vector<Bond::BondType>&>(),
             py::arg("atoms"), py::arg("bond_connectivity"), 
             py::arg("bond_types") = std::vector<Bond::BondType>(),
             "Create molecule with atoms and bonds")
        
        // Basic properties
        .def_property_readonly("atom_count", &Molecule::atom_count,
                              "Number of atoms in molecule")
        .def_property_readonly("bond_count", &Molecule::bond_count,
                              "Number of bonds in molecule")
        .def("is_empty", &Molecule::is_empty,
             "Check if molecule is empty (no atoms)")
        
        // Atom access and manipulation
        .def("atom", static_cast<const Atom& (Molecule::*)(std::size_t) const>(&Molecule::atom),
             py::arg("atom_index"), "Get atom by index",
             py::return_value_policy::reference)
        .def("atom", static_cast<Atom& (Molecule::*)(std::size_t)>(&Molecule::atom),
             py::arg("atom_index"), "Get mutable atom by index",
             py::return_value_policy::reference)
        .def("add_atom", static_cast<std::size_t (Molecule::*)(const Atom&)>(&Molecule::add_atom),
             py::arg("atom"), "Add atom to molecule (returns index)")
        .def("add_atom", static_cast<std::size_t (Molecule::*)(int, double, double, double)>(&Molecule::add_atom),
             py::arg("atomic_number"), py::arg("x"), py::arg("y"), py::arg("z"),
             "Add atom by atomic number and coordinates")
        .def("remove_atom", &Molecule::remove_atom, py::arg("atom_index"),
             "Remove atom and all associated bonds")
        
        // Bond access and manipulation
        .def("bond", static_cast<const Bond& (Molecule::*)(std::size_t) const>(&Molecule::bond),
             py::arg("bond_index"), "Get bond by index",
             py::return_value_policy::reference)
        .def("bond", static_cast<Bond& (Molecule::*)(std::size_t)>(&Molecule::bond),
             py::arg("bond_index"), "Get mutable bond by index",
             py::return_value_policy::reference)
        .def("add_bond", 
             static_cast<std::size_t (Molecule::*)(std::size_t, std::size_t, Bond::BondType)>(&Molecule::add_bond),
             py::arg("first_atom_index"), py::arg("second_atom_index"), 
             py::arg("bond_type") = Bond::BondType::SINGLE,
             "Add bond between atoms (returns bond index)")
        .def("add_bond",
             static_cast<std::size_t (Molecule::*)(std::size_t, std::size_t, double)>(&Molecule::add_bond),
             py::arg("first_atom_index"), py::arg("second_atom_index"), py::arg("bond_order"),
             "Add bond with explicit bond order")
        .def("remove_bond", &Molecule::remove_bond, py::arg("bond_index"),
             "Remove bond by index")
        .def("find_bond", &Molecule::find_bond,
             py::arg("first_atom_index"), py::arg("second_atom_index"),
             "Find bond connecting two atoms (returns index or None)")
        .def("bonds_for_atom", &Molecule::bonds_for_atom, py::arg("atom_index"),
             "Get all bond indices involving specified atom")
        
        // Molecular properties  
        .def("molecular_mass", &Molecule::molecular_mass,
             "Calculate total molecular mass (amu)")
        .def("center_of_mass", &Molecule::center_of_mass,
             "Calculate center of mass coordinates")
        .def("molecular_charge", &Molecule::molecular_charge,
             "Calculate total molecular charge")
        .def_property("molecular_multiplicity", &Molecule::molecular_multiplicity,
                     &Molecule::set_molecular_multiplicity,
                     "Spin multiplicity (2S + 1)")
        
        // Geometry operations
        .def("translate", &Molecule::translate, py::arg("displacement_vector"),
             "Translate entire molecule by vector")
        .def("center_at_origin", &Molecule::center_at_origin,
             "Translate molecule to place center of mass at origin")
        .def("bounding_box", &Molecule::bounding_box,
             "Get bounding box [min_x, min_y, min_z, max_x, max_y, max_z]")
        .def("molecular_diameter", &Molecule::molecular_diameter,
             "Calculate maximum distance between any two atoms")
        
        // Connectivity analysis
        .def("coordination_number", &Molecule::coordination_number,
             py::arg("atom_index"), "Get number of bonds for atom")
        .def("bonded_atoms", &Molecule::bonded_atoms, py::arg("atom_index"),
             "Get indices of atoms bonded to specified atom")
        .def("is_connected", &Molecule::is_connected,
             "Check if molecule is connected (single component)")
        .def("get_fragments", &Molecule::get_fragments,
             "Get separate molecular fragments as list of molecules")
        
        // Quantum chemistry compatibility
        .def("atomic_numbers", &Molecule::atomic_numbers,
             "Get atomic numbers as list")
        .def("element_symbols", &Molecule::element_symbols,
             "Get element symbols as list")
        .def("geometry_array", &Molecule::geometry_array,
             "Get flattened geometry array for QC software")
        .def("geometry_vector", &Molecule::geometry_vector,
             "Get geometry as list of coordinate arrays")
        .def("connectivity_matrix", &Molecule::connectivity_matrix,
             "Get bond order matrix")
        
        // File I/O and string representations
        .def("to_xyz_string", &Molecule::to_xyz_string,
             py::arg("comment") = "", "Generate XYZ format string")
        .def("to_detailed_string", &Molecule::to_detailed_string,
             "Generate detailed molecular information string")
        .def("to_string", &Molecule::to_string,
             "Brief string representation")
        .def("__str__", &Molecule::to_string)
        .def("__repr__", [](const Molecule& mol) {
             return "Molecule(" + std::to_string(mol.atom_count()) + " atoms, " +
                    std::to_string(mol.bond_count()) + " bonds)";
        })
        
        // Static construction methods
        .def_static("from_xyz_string", &Molecule::from_xyz_string,
                   py::arg("xyz_content"), "Create molecule from XYZ format string")
        .def_static("create_linear_molecule", &Molecule::create_linear_molecule,
                   py::arg("element_symbols"), py::arg("bond_length") = 1.5,
                   "Create linear molecule for testing")
        
        // Molecular property storage
        .def_property_readonly("properties", &Molecule::properties,
                              "Dictionary of molecular properties")
        .def("get_property", &Molecule::get_property, py::arg("property_name"),
             "Get molecular property by name")
        .def("set_property", &Molecule::set_property,
             py::arg("property_name"), py::arg("property_value"),
             "Set molecular property")
        
        // Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        
        // Python-specific convenience methods
        .def("__len__", &Molecule::atom_count, "Number of atoms (for len(molecule))")
        .def("__getitem__", [](const Molecule& mol, std::size_t index) -> const Atom& {
             return mol.atom(index);
        }, py::arg("index"), "Get atom by index (for molecule[i])",
           py::return_value_policy::reference)
        .def("__iter__", [](const Molecule& mol) {
             return py::make_iterator(
                 py::cast(&mol).attr("atom")(0),
                 py::cast(&mol).attr("atom")(mol.atom_count() - 1) + 1
             );
        }, py::keep_alive<0, 1>(), "Iterate over atoms");
    }

} // namespace libfrag
