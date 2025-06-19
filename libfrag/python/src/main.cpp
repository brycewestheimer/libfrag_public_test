#include "pybind11/pybind11.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"

#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <iostream>
#include <numeric>
#include <string>
#include <sstream>


// our headers
#include "libfrag/libfrag.hpp"
#include "libfrag/libfrag_config.hpp"

namespace py = pybind11;



namespace libfrag {


    // implementation in def_myclass.cpp
    void def_class(py::module & m);

    // implementation in def_build_config.cpp
    void def_build_config(py::module & m);

    // implementation in def_atom.cpp
    void def_atom(py::module & m);

    // implementation in def_bond.cpp
    void def_bond(py::module & m);

    // implementation in def_molecule.cpp
    void def_molecule(py::module & m);

    // implementation in def_fragment.cpp
    void def_fragment(py::module & m);

    // implementation in def_mbe.cpp
    void def_mbe(py::module & m);

}


// Python Module and Docstrings
PYBIND11_MODULE(_libfrag , module)
{
    xt::import_numpy();

    module.doc() = R"pbdoc(
        _libfrag  python bindings

        .. currentmodule:: _libfrag 

        .. autosummary::
           :toctree: _generate

           BuildConfiguration
           MyClass
           Atom
           Bond
           BondType
           Molecule
           Fragment
           FragmentLink
    )pbdoc";

    libfrag::def_build_config(module);
    libfrag::def_class(module);
    libfrag::def_atom(module);
    libfrag::def_bond(module);
    libfrag::def_molecule(module);
    libfrag::def_fragment(module);
    libfrag::def_mbe(module);

    // make version string
    std::stringstream ss;
    ss<<LIBFRAG_VERSION_MAJOR<<"."
      <<LIBFRAG_VERSION_MINOR<<"."
      <<LIBFRAG_VERSION_PATCH;
    module.attr("__version__") = ss.str().c_str();
}