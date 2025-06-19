#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/operators.h"
#include "pybind11/functional.h"

#include "xtensor-python/pyarray.hpp"

#include "libfrag/fragment.hpp"
#include "libfrag/bond.hpp"

namespace py = pybind11;

namespace libfrag {

    void def_fragment(py::module& m) {
        
        // First define FragmentLink class
        py::class_<FragmentLink>(m, "FragmentLink", R"pbdoc(
            Represents a covalent link between two fragments.
            
            FragmentLink stores information about bonds that connect different fragments,
            including the atoms involved and the bond properties.
            
            Examples
            --------
            >>> link = FragmentLink("frag1", "frag2", 0, 1, BondType.SINGLE, 1.0)
            >>> print(link.source_fragment_id)  # "frag1"
            >>> print(link.target_fragment_id)  # "frag2"
            >>> print(link.to_string())
        )pbdoc")
        
        .def(py::init<const std::string&, const std::string&, std::size_t, std::size_t, Bond::BondType, double>(),
             py::arg("source_fragment_id"), py::arg("target_fragment_id"), 
             py::arg("source_atom_index"), py::arg("target_atom_index"),
             py::arg("bond_type") = Bond::BondType::SINGLE, py::arg("bond_order") = 1.0,
             "Create a fragment link between two fragments")
        
        // Properties
        .def_property_readonly("source_fragment_id", &FragmentLink::source_fragment_id,
                              "ID of the source fragment")
        .def_property_readonly("target_fragment_id", &FragmentLink::target_fragment_id,
                              "ID of the target fragment")
        .def_property_readonly("source_atom_index", &FragmentLink::source_atom_index,
                              "Index of atom in source fragment")
        .def_property_readonly("target_atom_index", &FragmentLink::target_atom_index,
                              "Index of atom in target fragment")
        .def_property("bond_type", &FragmentLink::bond_type, &FragmentLink::set_bond_type,
                     "Type of bond connecting the fragments")
        .def_property("bond_order", &FragmentLink::bond_order, &FragmentLink::set_bond_order,
                     "Bond order (for fractional bonds)")
        
        // Methods
        .def("to_string", &FragmentLink::to_string, "String representation of the link")
        .def("__str__", &FragmentLink::to_string)
        .def("__repr__", [](const FragmentLink& link) {
            return "FragmentLink('" + link.source_fragment_id() + "' -> '" + 
                   link.target_fragment_id() + "')";
        })
        
        // Operators
        .def(py::self == py::self)
        .def(py::self != py::self);
          // Now define Fragment class - inherit from Molecule
        py::class_<Fragment, Molecule, std::shared_ptr<Fragment>>(m, "Fragment", R"pbdoc(
            Fragment class extending Molecule with fragment-specific functionality.
            
            The Fragment class implements the Composite design pattern, allowing fragments
            to contain subfragments. It extends Molecule to add fragment identification,
            inter-fragment connectivity, and hierarchical fragmentation capabilities.
            
            Key features:
            - Unique fragment identification system
            - Inter-fragment covalent links management
            - Composite pattern for hierarchical fragmentation
            - Fragment splitting and merging operations
            - Parent-child relationship tracking
            - Fragment-specific properties and metadata
            
            Examples
            --------
            >>> # Create atoms
            >>> atoms = [Atom(6, 0, 0, 0), Atom(1, 1, 0, 0), Atom(8, 3, 0, 0)]
            >>> bonds = [(0, 1), (1, 2)]
            >>> 
            >>> # Create fragment
            >>> fragment = Fragment(atoms, bonds, "test_fragment")
            >>> print(fragment.fragment_id)  # "test_fragment"
            >>> print(fragment.atom_count)   # 3
            >>> 
            >>> # Fragment by connectivity
            >>> subfragments = fragment.fragment_by_connectivity()
            >>> print(len(subfragments))
            >>> 
            >>> # Display fragment tree
            >>> print(fragment.to_tree_string())
        )pbdoc")
        
        // Inherit from Molecule (all Molecule methods will be available)
        // Note: pybind11 automatically handles inheritance
        
        // Constructors
        .def(py::init<>(), "Default constructor creates empty fragment")
        .def(py::init<const std::string&>(), py::arg("fragment_id"),
             "Construct fragment with unique ID")
        .def(py::init<const Molecule&, const std::string&>(),
             py::arg("molecule"), py::arg("fragment_id"),
             "Construct fragment from molecule with ID")
        .def(py::init<const std::vector<Atom>&, const std::string&>(),
             py::arg("atom_list"), py::arg("fragment_id"),
             "Construct fragment from atoms with ID")
        .def(py::init<const std::vector<Atom>&, 
                      const std::vector<std::pair<Fragment::AtomIndex, Fragment::AtomIndex>>&,
                      const std::string&,
                      const std::vector<Bond::BondType>&>(),
             py::arg("atom_list"), py::arg("bond_connectivity"), py::arg("fragment_id"),
             py::arg("bond_types") = std::vector<Bond::BondType>(),
             "Construct fragment with atoms, bonds, and ID")
        
        // Fragment identification
        .def_property("fragment_id", &Fragment::fragment_id, &Fragment::set_fragment_id,
                     "Unique fragment identifier")
        .def_property_readonly("generation_level", &Fragment::generation_level,
                              "Generation level in fragment hierarchy")
        
        // Parent-child relationship management
        .def_property_readonly("parent_fragment", 
                              [](const Fragment& f) -> py::object {
                                  auto parent = f.parent_fragment().lock();
                                  if (parent) {
                                      return py::cast(parent);
                                  } else {
                                      return py::none();
                                  }
                              },
                              "Parent fragment (None if root)")
        .def("set_parent_fragment", &Fragment::set_parent_fragment, py::arg("parent"),
             "Set parent fragment")
        .def("is_root_fragment", &Fragment::is_root_fragment,
             "Check if this is a root fragment (no parent)")
        .def("get_root_fragment", &Fragment::get_root_fragment,
             "Get root fragment in hierarchy")
        
        // Subfragment management (Composite pattern)
        .def_property_readonly("subfragment_count", &Fragment::subfragment_count,
                              "Number of direct subfragments")
        .def_property_readonly("subfragments", &Fragment::subfragments,
                              "List of direct subfragments")
        .def("subfragment", &Fragment::subfragment, py::arg("index"),
             "Get subfragment by index")
        .def("add_subfragment", &Fragment::add_subfragment, py::arg("subfragment"),
             "Add subfragment to this fragment")
        .def("remove_subfragment", &Fragment::remove_subfragment, py::arg("index"),
             "Remove subfragment by index")
        .def("clear_subfragments", &Fragment::clear_subfragments,
             "Remove all subfragments")
        .def("has_subfragments", &Fragment::has_subfragments,
             "Check if fragment has subfragments")
        .def("get_all_subfragments", &Fragment::get_all_subfragments,
             "Get all fragments in subtree (recursive)")
        .def("total_atom_count", &Fragment::total_atom_count,
             "Get total number of atoms in entire fragment tree")
        
        // Inter-fragment connectivity
        .def_property_readonly("fragment_links", &Fragment::fragment_links,
                              "List of fragment links")
        .def("add_fragment_link", &Fragment::add_fragment_link, py::arg("link"),
             "Add link to another fragment")
        .def("remove_fragment_link", &Fragment::remove_fragment_link, py::arg("index"),
             "Remove fragment link by index")
        .def("find_links_to_fragment", &Fragment::find_links_to_fragment, 
             py::arg("target_fragment_id"), "Find links to specific fragment")
        .def("linked_fragment_ids", &Fragment::linked_fragment_ids,
             "Get all fragment IDs this fragment is linked to")
        .def("is_linked_to", &Fragment::is_linked_to, py::arg("target_fragment_id"),
             "Check if fragment is linked to another fragment")
        
        // Fragment splitting operations
        .def("fragment_by_connectivity", &Fragment::fragment_by_connectivity,
             py::arg("max_fragment_size") = 0,
             "Split fragment into subfragments based on connectivity")
        .def("fragment_by_bond_breaking", &Fragment::fragment_by_bond_breaking,
             py::arg("bonds_to_break"),
             "Split fragment by breaking specific bonds")
        .def("fragment_by_custom_function", &Fragment::fragment_by_custom_function,
             py::arg("fragmenter"),
             "Split fragment using custom fragmentation function")
        .def("fragment_by_size", &Fragment::fragment_by_size,
             py::arg("target_fragment_count"),
             "Split fragment into roughly equal-sized pieces")
        
        // Fragment merging operations
        .def("merge_subfragments", &Fragment::merge_subfragments,
             py::arg("preserve_links") = true,
             "Merge all subfragments back into this fragment")
        .def("merge_selected_subfragments", &Fragment::merge_selected_subfragments,
             py::arg("subfragment_indices"), py::arg("preserve_links") = true,
             "Merge specific subfragments")
        
        // Fragment analysis
        .def("calculate_fragment_properties", &Fragment::calculate_fragment_properties,
             "Calculate fragment-specific properties")
        .def("fragment_composition", &Fragment::fragment_composition,
             "Get fragment composition (element counts)")
        .def("fragment_hash", &Fragment::fragment_hash,
             "Generate fragment hash for quick comparison")
        
        // String representations
        .def("to_fragment_string", &Fragment::to_fragment_string,
             "Generate detailed fragment information")
        .def("to_tree_string", &Fragment::to_tree_string,
             py::arg("indent_level") = 0,
             "Generate fragment tree representation")
        .def("__str__", &Fragment::to_fragment_string)
        .def("__repr__", [](const Fragment& f) {
            return "Fragment(id='" + f.fragment_id() + "', atoms=" + 
                   std::to_string(f.atom_count()) + ", subfragments=" + 
                   std::to_string(f.subfragment_count()) + ")";
        })
        
        // Static factory methods
        .def_static("create_from_molecule", &Fragment::create_from_molecule,
                   py::arg("molecule"), py::arg("id_prefix") = "frag",
                   "Create fragment from molecule with automatic ID generation")
        .def_static("create_fragment_hierarchy", &Fragment::create_fragment_hierarchy,
                   py::arg("molecule"), py::arg("max_fragment_size") = 10,
                   py::arg("id_prefix") = "frag",
                   "Create fragment hierarchy from molecule using automatic fragmentation")
        
        // Operators
        .def(py::self == py::self)
        .def(py::self != py::self);
        
        // Add some useful utility functions
        m.def("create_test_fragment", []() {
            std::vector<Atom> atoms = {
                Atom(6, 0.0, 0.0, 0.0),    // Carbon
                Atom(1, 1.1, 0.0, 0.0),    // Hydrogen
                Atom(1, -0.5, 0.9, 0.0),   // Hydrogen
                Atom(1, -0.5, -0.9, 0.0),  // Hydrogen
            };
            
            std::vector<std::pair<Fragment::AtomIndex, Fragment::AtomIndex>> bonds = {
                {0, 1}, {0, 2}, {0, 3}
            };
            
            return std::make_shared<Fragment>(atoms, bonds, "test_methane");
        }, "Create a test methane fragment for demonstration");
        
        m.def("create_complex_test_fragment", []() {
            std::vector<Atom> atoms = {
                // Methane part
                Atom(6, 0.0, 0.0, 0.0),    // Carbon
                Atom(1, 1.1, 0.0, 0.0),    // H1
                Atom(1, -0.5, 0.9, 0.0),   // H2
                Atom(1, -0.5, -0.9, 0.0),  // H3
                
                // Water part (disconnected from methane)
                Atom(8, 5.0, 0.0, 0.0),    // Oxygen
                Atom(1, 5.5, 0.8, 0.0),    // H4
                Atom(1, 5.5, -0.8, 0.0),   // H5
            };
            
            std::vector<std::pair<Fragment::AtomIndex, Fragment::AtomIndex>> bonds = {
                {0, 1}, {0, 2}, {0, 3},  // Methane bonds
                {4, 5}, {4, 6}           // Water bonds
            };
            
            return std::make_shared<Fragment>(atoms, bonds, "test_complex");
        }, "Create a complex test fragment with disconnected components");
    }

} // namespace libfrag
