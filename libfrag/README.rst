=================================================
libfrag
=================================================

.. image:: https://readthedocs.org/projects/libfrag/badge/?version=latest
        :target: http://libfrag.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://github.com/yourusername/libfrag/workflows/CI/badge.svg
        :target: https://github.com/yourusername/libfrag/actions

.. image:: https://img.shields.io/travis/brycewestheimer/libfrag.svg
        :target: https://travis-ci.org/brycewestheimer/libfrag

.. image:: https://travis-ci.org/brycewestheimer/libfrag.svg?branch=master
    :target: https://travis-ci.org/brycewestheimer/libfrag

.. image:: https://circleci.com/gh/brycewestheimer/libfrag/tree/master.svg?style=svg
    :target: https://circleci.com/gh/brycewestheimer/libfrag/tree/master

.. image:: https://dev.azure.com/brycewestheimer/libfrag/_apis/build/status/brycewestheimer.libfrag?branchName=master
    :target: https://dev.azure.com/brycewestheimer/libfrag/_build/latest?definitionId=1&branchName=master


Overview
--------

libfrag is a modern C++ library for quantum chemistry calculations using fragment-based methods. It provides both C++ and Python APIs for performing efficient molecular fragmentation and analysis.

Features
--------

Current features include:

* Modern C++17 implementation
* Build system with CMake 3.15+
* Python bindings using pybind11
* Integration with xtensor for efficient numerical computations
* Cross-platform support (Linux, macOS, Windows)
* Comprehensive test suite
* Documentation with Sphinx

Installation
------------

From source::

    git clone https://github.com/yourusername/libfrag.git
    cd libfrag
    mkdir build && cd build
    cmake ..
    make -j4

Python Package::

    pip install libfrag

Quick Start
-----------

C++::

    #include "libfrag/libfrag.hpp"
    
    int main() {
        libfrag::MyClass frag(42);
        frag.hello_world();
        return 0;
    }

Python::

    import libfrag
    
    # Use the library
    result = libfrag.pure_python()
    print(result)

License
-------

This project is licensed under the MIT License - see the LICENSE.txt file for details.