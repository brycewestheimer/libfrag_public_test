#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/operators.h"

#include "xtensor-python/pyarray.hpp"

#include "libfrag/atom.hpp"

namespace py = pybind11;

namespace libfrag {

    void def_atom(py::module& m) {
        py::class_<Atom>(m, "Atom", R"pbdoc(
            Comprehensive Atom class for quantum chemistry applications.
            
            Designed to interface with MolSSI Driver Interface, QCEngine, and QCElemental.
            Stores atomic properties including coordinates, charge, mass, and quantum properties.
            
            Examples
            --------
            Create atoms in different ways:
            
            >>> # From atomic number and coordinates
            >>> atom1 = Atom(6, 0.0, 0.0, 0.0)  # Carbon at origin
            
            >>> # From element symbol and coordinates  
            >>> atom2 = Atom("N", 1.0, 0.0, 0.0)  # Nitrogen
            
            >>> # From numpy array
            >>> import numpy as np
            >>> coords = np.array([0.0, 1.0, 0.0])
            >>> atom3 = Atom(8, coords)  # Oxygen
            
            >>> # Access properties
            >>> print(atom1.symbol)  # "C"
            >>> print(atom1.atomic_mass)  # 12.011
            >>> print(atom1.distance_to(atom2))  # 1.0
        )pbdoc")
        
        // Constructors
        .def(py::init<>(), "Default constructor (creates hydrogen at origin)")
        .def(py::init<int, double, double, double>(), 
             py::arg("atomic_number"), py::arg("x"), py::arg("y"), py::arg("z"),
             "Create atom from atomic number and coordinates")
        .def(py::init<const std::string&, double, double, double>(),
             py::arg("symbol"), py::arg("x"), py::arg("y"), py::arg("z"), 
             "Create atom from element symbol and coordinates")
        .def(py::init<int, const std::array<double, 3>&>(),
             py::arg("atomic_number"), py::arg("coordinates"),
             "Create atom from atomic number and coordinate array")
        .def(py::init<int, const xt::xarray<double>&>(),
             py::arg("atomic_number"), py::arg("coordinates"),
             "Create atom from atomic number and numpy array")
        
        // Properties (read-only)
        .def_property_readonly("atomic_number", &Atom::atomic_number,
                              "Atomic number (Z)")
        .def_property_readonly("symbol", &Atom::symbol,
                              "Element symbol (e.g., 'C', 'N', 'O')")
        .def_property_readonly("coordinates", &Atom::coordinates,
                              "Coordinates as [x, y, z] array")
        .def_property_readonly("x", &Atom::x, "X coordinate")
        .def_property_readonly("y", &Atom::y, "Y coordinate") 
        .def_property_readonly("z", &Atom::z, "Z coordinate")
        .def_property_readonly("atomic_mass", &Atom::atomic_mass,
                              "Atomic mass in atomic mass units")
        .def_property_readonly("mass_number", &Atom::mass_number,
                              "Mass number of isotope")
        
        // Charge and electronic properties
        .def_property("formal_charge", &Atom::formal_charge, &Atom::set_formal_charge,
                     "Formal charge")
        .def_property("partial_charge", &Atom::partial_charge, &Atom::set_partial_charge,
                     "Partial charge (e.g., from population analysis)")
        .def_property("multiplicity", &Atom::multiplicity, &Atom::set_multiplicity,
                     "Spin multiplicity")
        
        // Quantum properties
        .def_property_readonly("properties", &Atom::properties,
                              "Dictionary of quantum mechanical properties")
        
        // Methods for coordinates
        .def("set_coordinates", 
             static_cast<void (Atom::*)(double, double, double)>(&Atom::set_coordinates),
             py::arg("x"), py::arg("y"), py::arg("z"),
             "Set coordinates from x, y, z values")
        .def("set_coordinates",
             static_cast<void (Atom::*)(const std::array<double, 3>&)>(&Atom::set_coordinates),
             py::arg("coordinates"),
             "Set coordinates from array")
        .def("set_coordinates",
             static_cast<void (Atom::*)(const xt::xarray<double>&)>(&Atom::set_coordinates),
             py::arg("coordinates"),
             "Set coordinates from numpy array")
        
        // Property management
        .def("get_property", &Atom::get_property, py::arg("key"),
             "Get a quantum property by key (returns None if not found)")
        .def("set_property", &Atom::set_property, py::arg("key"), py::arg("value"),
             "Set a quantum property")
        .def("set_mass_number", &Atom::set_mass_number, py::arg("mass_number"),
             "Set the isotope mass number")
        
        // Utility methods
        .def("distance_to", &Atom::distance_to, py::arg("other"),
             "Calculate distance to another atom")
        .def("vector_to", &Atom::vector_to, py::arg("other"),
             "Calculate vector from this atom to another")
        
        // Conversion methods
        .def("coordinates_array", &Atom::coordinates_array,
             "Get coordinates as numpy array")
        .def("coordinates_vector", &Atom::coordinates_vector,
             "Get coordinates as list")
        
        // String representations
        .def("to_string", &Atom::to_string,
             "String representation of atom")
        .def("to_xyz_string", &Atom::to_xyz_string,
             "XYZ format string representation")
        .def("__str__", &Atom::to_string)
        .def("__repr__", [](const Atom& atom) {
            return "Atom(" + std::to_string(atom.atomic_number()) + ", " +
                   std::to_string(atom.x()) + ", " + std::to_string(atom.y()) + ", " +
                   std::to_string(atom.z()) + ")";
        })
        
        // Operators
        .def(py::self == py::self)
        .def(py::self != py::self)
        
        // Static methods
        .def_static("symbol_to_atomic_number", &Atom::symbol_to_atomic_number,
                   py::arg("symbol"), "Convert element symbol to atomic number")
        .def_static("atomic_number_to_symbol", &Atom::atomic_number_to_symbol,
                   py::arg("atomic_number"), "Convert atomic number to element symbol")
        .def_static("is_valid_atomic_number", &Atom::is_valid_atomic_number,
                   py::arg("atomic_number"), "Check if atomic number is valid")
        .def_static("is_valid_element_symbol", &Atom::is_valid_element_symbol,
                   py::arg("symbol"), "Check if element symbol is valid");
        
        // Add some useful class attributes for QC compatibility
        m.attr("BOHR_TO_ANGSTROM") = 0.52917721067;
        m.attr("ANGSTROM_TO_BOHR") = 1.8897261246;
        m.attr("HARTREE_TO_EV") = 27.211386245988;
        m.attr("EV_TO_HARTREE") = 0.036749322176;
    }

} // namespace libfrag
