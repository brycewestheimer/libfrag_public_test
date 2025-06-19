#!/usr/bin/env python3
"""
Many-Body Expansion (MBE) Example for libfrag

This example demonstrates how to use the MBE functionality in libfrag
to perform fragment-based quantum chemistry calculations.
"""

import libfrag
from libfrag.mbe import *
import numpy as np

def main():
    print("=== libfrag Many-Body Expansion Demo ===")
    
    # Create a simple water dimer system
    atoms = [
        libfrag.Atom("O", 0.0, 0.0, 0.0),      # Water 1
        libfrag.Atom("H", 0.76, 0.59, 0.0),
        libfrag.Atom("H", -0.76, 0.59, 0.0),
        libfrag.Atom("O", 3.0, 0.0, 0.0),      # Water 2  
        libfrag.Atom("H", 3.76, 0.59, 0.0),
        libfrag.Atom("H", 2.24, 0.59, 0.0)
    ]
    
    molecule = libfrag.Molecule(atoms)
    print(f"Created water dimer with {len(molecule)} atoms")
    
    # Configure MBE calculation
    config = MBEConfig.default_2body()
    config.set_qm_method("HF")
    config.set_basis_set("6-31G*")
    config.set_fragmentation_scheme(FragmentationScheme.MOLECULAR)
    config.set_distance_cutoff(8.0)
    
    print("\nMBE Configuration:")
    print(config.to_string())
    
    # Test fragment generation
    print("\n=== Fragment Generation ===")
    generator = MBEFragmentGenerator(config)
    
    # Estimate number of combinations
    n_fragments_estimate = 2  # Two water molecules
    for order in range(1, config.max_order() + 1):
        n_combinations = generator.estimate_combinations(n_fragments_estimate, order)
        print(f"{order}-body combinations: {n_combinations}")
    
    # Test MBE results storage
    print("\n=== Results Testing ===")
    results = MBEResults(config)
    
    # Add some mock fragment results
    for i in range(2):  # 1-body terms
        fragment_result = FragmentCalculationResult([i], 1)
        fragment_result.total_energy = -76.0 + i * 0.001  # Mock energies
        fragment_result.qm_method = "HF"
        fragment_result.basis_set = "6-31G*"
        fragment_result.converged = True
        results.add_fragment_result(fragment_result)
    
    # Add 2-body interaction
    fragment_result = FragmentCalculationResult([0, 1], 2)
    fragment_result.total_energy = -152.005  # Mock interaction energy
    fragment_result.qm_method = "HF"
    fragment_result.basis_set = "6-31G*"
    fragment_result.converged = True
    results.add_fragment_result(fragment_result)
    
    print(f"Added {results.n_calculations()} fragment calculations")
    print(f"Total MBE energy: {results.total_energy():.6f} Hartree")
    
    # Print energy decomposition
    print("\nEnergy Contributions:")
    for order in range(1, results.max_order() + 1):
        contribution = results.energy_contribution(order)
        print(f"  {order}-body: {contribution:.6f} Hartree")
    
    # Export results
    print("\n=== Results Export ===")
    print("Summary Report:")
    print(results.summary_report())
    
    print("\nJSON Export (first 200 chars):")
    json_output = results.to_json()
    print(json_output[:200] + "..." if len(json_output) > 200 else json_output)
    
    # Test different configurations
    print("\n=== Configuration Examples ===")
    
    # 3-body configuration
    config_3body = MBEConfig.default_3body()
    print(f"3-body config: max_order={config_3body.max_order()}")
    
    # Custom fragmentation
    custom_config = MBEConfig(2, FragmentationScheme.CUSTOM)
    custom_fragments = [
        [0, 1, 2],  # First water
        [3, 4, 5]   # Second water  
    ]
    custom_config.set_custom_fragments(custom_fragments)
    custom_config.set_qm_method("B3LYP")
    custom_config.set_basis_set("def2-TZVP")
    
    print(f"Custom config: {len(custom_config.custom_fragments())} fragments defined")
    print(f"QM Method: {custom_config.qm_method()}/{custom_config.basis_set()}")
    
    # Test utility functions
    print("\n=== Utility Functions ===")
    n_combinations_2of5 = binomial_coefficient(5, 2)
    print(f"Binomial coefficient C(5,2) = {n_combinations_2of5}")
    
    combinations_3of4 = generate_combinations(4, 3)
    print(f"All 3-combinations from 4 items: {combinations_3of4}")
    
    # Test overlap checking
    combo1 = [0, 1, 2]
    combo2 = [2, 3, 4]
    overlap = combinations_overlap(combo1, combo2)
    print(f"Combinations {combo1} and {combo2} overlap: {overlap}")
    
    print("\n=== Demo Complete ===")
    print("MBE functionality successfully demonstrated!")
    print("\nKey features showcased:")
    print("- Configuration management with different fragmentation schemes")
    print("- Results storage and energy decomposition analysis")
    print("- Fragment combination generation and utilities")
    print("- Export capabilities (JSON, CSV, summary reports)")
    print("- Validation and error checking")
    
    print("\nNote: This demo uses mock data since QM calculator integration")
    print("requires connection to actual quantum chemistry software.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        print("Note: This example requires the libfrag Python module to be built and installed.")
        print("Run this after successful compilation of the C++ library with Python bindings.")
