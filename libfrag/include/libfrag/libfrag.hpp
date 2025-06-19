#pragma once
#ifndef LIBFRAG_LIBFRAG_HPP
#define LIBFRAG_LIBFRAG_HPP

#include <cstdint>
#include <iostream>

// Include all libfrag components
#include "libfrag/atom.hpp"
#include "libfrag/bond.hpp"
#include "libfrag/molecule.hpp"
#include "libfrag/fragment.hpp"

// Include Many-Body Expansion components
#include "libfrag/mbe_config.hpp"
#include "libfrag/mbe_results.hpp"
#include "libfrag/mbe_fragment_generator.hpp"
#include "libfrag/mbe_calculator.hpp"

// Include Global SCF components
#include "libfrag/global_scf_config.hpp"
#include "libfrag/global_scf_results.hpp"
#include "libfrag/global_scf_fragment_generator.hpp"
#include "libfrag/global_scf_calculator.hpp"

namespace libfrag {
    
    class MyClass
    {
    public:
        MyClass(const uint64_t size)
        : m_size(size)
        {

        }
        
        void hello_world()
        {
            std::cout<<"Hello World!\n";
        }
    private:
        uint64_t m_size;
    };

} // end namespace libfrag


#endif // LIBFRAG_LIBFRAG_HPP