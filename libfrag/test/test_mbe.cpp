#include "libfrag/mbe_config.hpp"
#include "libfrag/mbe_results.hpp"
#include "libfrag/mbe_fragment_generator.hpp"
#include "libfrag/mbe_calculator.hpp"
#include <iostream>

int main() {
    using namespace libfrag;
    
    std::cout << "Testing MBE header compilation..." << std::endl;
    
    // Test MBEConfig
    MBEConfig config(2, MBEConfig::FragmentationScheme::MOLECULAR);
    config.set_qm_method("HF");
    config.set_basis_set("6-31G*");
    std::cout << "MBEConfig created successfully" << std::endl;
    
    // Test MBEResults
    MBEResults results;
    FragmentCalculationResult fragment_result;
    fragment_result.fragment_id = "test_fragment";
    fragment_result.n_body_order = 1;
    fragment_result.total_energy = -1.0;
    fragment_result.converged = true;
    results.add_fragment_result(fragment_result);
    std::cout << "MBEResults created successfully" << std::endl;
    
    // Test MBEFragmentGenerator
    MBEFragmentGenerator generator(config);
    std::cout << "MBEFragmentGenerator created successfully" << std::endl;
    
    std::cout << "All MBE headers compiled successfully!" << std::endl;
    
    return 0;
}
