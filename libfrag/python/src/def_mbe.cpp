#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/chrono.h"
#include "pybind11/functional.h"
#include "xtensor-python/pyarray.hpp"

#include "libfrag/mbe_config.hpp"
#include "libfrag/mbe_results.hpp"
#include "libfrag/mbe_fragment_generator.hpp"
#include "libfrag/mbe_calculator.hpp"

namespace py = pybind11;
using namespace libfrag;

namespace libfrag {

    void def_mbe_config(py::module & m) {
        py::enum_<MBEConfig::FragmentationScheme>(m, "FragmentationScheme")
            .value("ATOMIC", MBEConfig::FragmentationScheme::ATOMIC)
            .value("MOLECULAR", MBEConfig::FragmentationScheme::MOLECULAR)
            .value("FUNCTIONAL_GROUP", MBEConfig::FragmentationScheme::FUNCTIONAL_GROUP)
            .value("CUSTOM", MBEConfig::FragmentationScheme::CUSTOM)
            .value("DISTANCE_BASED", MBEConfig::FragmentationScheme::DISTANCE_BASED);
            
        py::enum_<MBEConfig::TruncationMethod>(m, "TruncationMethod")
            .value("ORDER_BASED", MBEConfig::TruncationMethod::ORDER_BASED)
            .value("DISTANCE_BASED", MBEConfig::TruncationMethod::DISTANCE_BASED)
            .value("ENERGY_BASED", MBEConfig::TruncationMethod::ENERGY_BASED);

        py::class_<MBEConfig>(m, "MBEConfig")
            .def(py::init<>())
            .def(py::init<int, MBEConfig::FragmentationScheme>(),
                 py::arg("max_order"), py::arg("scheme") = MBEConfig::FragmentationScheme::MOLECULAR)
            
            // Configuration setters
            .def("set_max_order", &MBEConfig::set_max_order, py::arg("order"))
            .def("set_fragmentation_scheme", &MBEConfig::set_fragmentation_scheme, py::arg("scheme"))
            .def("set_distance_cutoff", &MBEConfig::set_distance_cutoff, py::arg("cutoff"))
            .def("set_energy_threshold", &MBEConfig::set_energy_threshold, py::arg("threshold"))
            .def("set_qm_method", &MBEConfig::set_qm_method, py::arg("method"))
            .def("set_basis_set", &MBEConfig::set_basis_set, py::arg("basis"))
            .def("set_charge_embedding", &MBEConfig::set_charge_embedding, py::arg("enable"))
            .def("set_custom_fragments", &MBEConfig::set_custom_fragments, py::arg("fragments"))
            .def("set_qm_options", &MBEConfig::set_qm_options, py::arg("options"))
            
            // Configuration getters
            .def("max_order", &MBEConfig::max_order)
            .def("fragmentation_scheme", &MBEConfig::fragmentation_scheme)
            .def("distance_cutoff", &MBEConfig::distance_cutoff)
            .def("energy_threshold", &MBEConfig::energy_threshold)
            .def("qm_method", &MBEConfig::qm_method)
            .def("basis_set", &MBEConfig::basis_set)
            .def("charge_embedding", &MBEConfig::charge_embedding)
            .def("truncation_method", &MBEConfig::truncation_method)
            .def("custom_fragments", &MBEConfig::custom_fragments)
            .def("qm_options", &MBEConfig::qm_options)
            
            // Utilities
            .def("is_valid", &MBEConfig::is_valid)
            .def("to_string", &MBEConfig::to_string)
            .def_static("from_string", &MBEConfig::from_string, py::arg("config_string"))
            .def_static("default_2body", &MBEConfig::default_2body)
            .def_static("default_3body", &MBEConfig::default_3body)
            
            .def("__repr__", [](const MBEConfig& config) {
                return "<MBEConfig max_order=" + std::to_string(config.max_order()) + 
                       " method=" + config.qm_method() + "/" + config.basis_set() + ">";
            });
    }

    void def_mbe_results(py::module & m) {
        py::class_<FragmentCalculationResult>(m, "FragmentCalculationResult")
            .def(py::init<>())
            .def(py::init<const std::vector<std::size_t>&, int>(),
                 py::arg("fragment_indices"), py::arg("n_body_order"))
            
            // Data members
            .def_readwrite("fragment_indices", &FragmentCalculationResult::fragment_indices)
            .def_readwrite("n_body_order", &FragmentCalculationResult::n_body_order)
            .def_readwrite("fragment_id", &FragmentCalculationResult::fragment_id)
            .def_readwrite("total_energy", &FragmentCalculationResult::total_energy)
            .def_readwrite("nuclear_repulsion", &FragmentCalculationResult::nuclear_repulsion)
            .def_readwrite("electronic_energy", &FragmentCalculationResult::electronic_energy)
            .def_readwrite("correlation_energy", &FragmentCalculationResult::correlation_energy)
            .def_readwrite("properties", &FragmentCalculationResult::properties)
            .def_readwrite("qm_method", &FragmentCalculationResult::qm_method)
            .def_readwrite("basis_set", &FragmentCalculationResult::basis_set)
            .def_readwrite("n_electrons", &FragmentCalculationResult::n_electrons)
            .def_readwrite("multiplicity", &FragmentCalculationResult::multiplicity)
            .def_readwrite("charge", &FragmentCalculationResult::charge)
            .def_readwrite("computation_time", &FragmentCalculationResult::computation_time)
            .def_readwrite("converged", &FragmentCalculationResult::converged)
            .def_readwrite("scf_iterations", &FragmentCalculationResult::scf_iterations)
            .def_readwrite("status", &FragmentCalculationResult::status)
            
            // Methods
            .def("to_string", &FragmentCalculationResult::to_string)
            .def("is_valid", &FragmentCalculationResult::is_valid)
            
            .def("__repr__", [](const FragmentCalculationResult& result) {
                return "<FragmentCalculationResult " + result.fragment_id + 
                       " energy=" + std::to_string(result.total_energy) + ">";
            });

        py::class_<MBEResults>(m, "MBEResults")
            .def(py::init<>())
            .def(py::init<const MBEConfig&>(), py::arg("config"))
            
            // Result management
            .def("add_fragment_result", 
                 py::overload_cast<const FragmentCalculationResult&>(&MBEResults::add_fragment_result),
                 py::arg("result"))
            
            // Energy accessors
            .def("total_energy", &MBEResults::total_energy)
            .def("energy_contribution", &MBEResults::energy_contribution, py::arg("order"))
            .def("energy_contributions", &MBEResults::energy_contributions)
            .def("fragment_energies", &MBEResults::fragment_energies)
            .def("interaction_energies", &MBEResults::interaction_energies, py::arg("order"))
            
            // Result accessors
            .def("fragment_results", &MBEResults::fragment_results)
            .def("results_by_order", &MBEResults::results_by_order, py::arg("order"))
            .def("max_order", &MBEResults::max_order)
            .def("n_fragments", &MBEResults::n_fragments)
            
            // Error analysis
            .def("estimated_truncation_error", &MBEResults::estimated_truncation_error)
            .def("is_converged", &MBEResults::is_converged, 
                 py::arg("threshold") = 1e-6)
            .def("convergence_analysis", &MBEResults::convergence_analysis)
            
            // Timing
            .def("total_computation_time", &MBEResults::total_computation_time)
            .def("computation_time_by_order", &MBEResults::computation_time_by_order, py::arg("order"))
            .def("performance_statistics", &MBEResults::performance_statistics)
            
            // Export
            .def("to_json", &MBEResults::to_json)
            .def("to_csv", &MBEResults::to_csv)
            .def("summary_report", &MBEResults::summary_report)
            .def("energy_decomposition_table", &MBEResults::energy_decomposition_table)
            
            // Utilities
            .def("validate", &MBEResults::validate)
            .def("clear", &MBEResults::clear)
            .def("empty", &MBEResults::empty)
            .def("n_calculations", &MBEResults::n_calculations)
            
            .def("__len__", &MBEResults::n_calculations)
            .def("__getitem__", &MBEResults::operator[], py::arg("index"))
            .def("__repr__", [](const MBEResults& results) {
                return "<MBEResults total_energy=" + std::to_string(results.total_energy()) + 
                       " max_order=" + std::to_string(results.max_order()) + 
                       " n_calculations=" + std::to_string(results.n_calculations()) + ">";
            });
    }

    void def_mbe_fragment_generator(py::module & m) {
        py::class_<MBEFragmentGenerator>(m, "MBEFragmentGenerator")
            .def(py::init<>())
            .def(py::init<const MBEConfig&>(), py::arg("config"))
            
            // Main generation methods
            .def("generate_all_combinations", &MBEFragmentGenerator::generate_all_combinations,
                 py::arg("molecule"), py::arg("config"))
            .def("generate_n_body_combinations", &MBEFragmentGenerator::generate_n_body_combinations,
                 py::arg("fragments"), py::arg("n_body"), py::arg("config"))
            
            // Fragmentation methods
            .def("create_initial_fragments", &MBEFragmentGenerator::create_initial_fragments,
                 py::arg("molecule"), py::arg("config"))
            .def("fragment_by_atoms", &MBEFragmentGenerator::fragment_by_atoms, py::arg("molecule"))
            .def("fragment_by_molecules", &MBEFragmentGenerator::fragment_by_molecules, py::arg("molecule"))
            .def("fragment_by_functional_groups", &MBEFragmentGenerator::fragment_by_functional_groups, py::arg("molecule"))
            .def("fragment_by_distance", &MBEFragmentGenerator::fragment_by_distance,
                 py::arg("molecule"), py::arg("cutoff"))
            .def("fragment_by_custom", &MBEFragmentGenerator::fragment_by_custom,
                 py::arg("molecule"), py::arg("fragment_definitions"))
            
            // Filtering
            .def("filter_by_distance", &MBEFragmentGenerator::filter_by_distance,
                 py::arg("combinations"), py::arg("fragments"), py::arg("cutoff"))
            .def("filter_by_energy", &MBEFragmentGenerator::filter_by_energy,
                 py::arg("combinations"), py::arg("fragments"), py::arg("threshold"))
            .def("remove_duplicates", &MBEFragmentGenerator::remove_duplicates, py::arg("combinations"))
            
            // Utilities
            .def("calculate_fragment_distance", &MBEFragmentGenerator::calculate_fragment_distance,
                 py::arg("frag1"), py::arg("frag2"))
            .def("is_valid_combination", &MBEFragmentGenerator::is_valid_combination,
                 py::arg("combination"), py::arg("fragments"), py::arg("config"))
            .def("estimate_combinations", &MBEFragmentGenerator::estimate_combinations,
                 py::arg("n_fragments"), py::arg("n_body"))
            .def("generation_statistics", &MBEFragmentGenerator::generation_statistics)
            
            // Configuration
            .def("set_config", &MBEFragmentGenerator::set_config, py::arg("config"))
            .def("config", &MBEFragmentGenerator::config)
            
            .def("__repr__", [](const MBEFragmentGenerator& gen) {
                return "<MBEFragmentGenerator>";
            });
        
        // Utility functions
        m.def("binomial_coefficient", &mbe_utils::binomial_coefficient,
              py::arg("n"), py::arg("k"));
        m.def("generate_combinations", &mbe_utils::generate_combinations,
              py::arg("n"), py::arg("k"));
        m.def("combinations_overlap", &mbe_utils::combinations_overlap,
              py::arg("combo1"), py::arg("combo2"));
    }

    void def_mbe_calculator(py::module & m) {
        // QM Calculator Interface (abstract base class)
        py::class_<QMCalculatorInterface>(m, "QMCalculatorInterface")
            .def("calculate_energy", &QMCalculatorInterface::calculate_energy,
                 py::arg("fragment"), py::arg("method"), py::arg("basis_set"),
                 py::arg("charge") = 0, py::arg("multiplicity") = 1,
                 py::arg("options") = std::unordered_map<std::string, std::string>{})
            .def("is_available", &QMCalculatorInterface::is_available)
            .def("software_info", &QMCalculatorInterface::software_info);

        py::class_<MBECalculator>(m, "MBECalculator")
            .def(py::init<>())
            // Note: Constructor with QMCalculatorPtr requires special handling for Python
            .def(py::init<std::unique_ptr<QMCalculatorInterface>, const MBEConfig&>(),
                 py::arg("qm_calculator"), py::arg("config"))
            
            // Configuration
            .def("set_config", &MBECalculator::set_config, py::arg("config"))
            .def("config", &MBECalculator::config)
            
            // Callbacks (simplified for Python)
            .def("set_progress_callback", 
                 [](MBECalculator& calc, py::function callback) {
                     calc.set_progress_callback([callback](double progress, const std::string& status) {
                         callback(progress, status);
                     });
                 }, py::arg("callback"))
            .def("set_log_callback",
                 [](MBECalculator& calc, py::function callback) {
                     calc.set_log_callback([callback](const std::string& message) {
                         callback(message);
                     });
                 }, py::arg("callback"))
            
            // Main calculation methods
            .def("calculate", 
                 py::overload_cast<const Molecule&>(&MBECalculator::calculate),
                 py::arg("molecule"))
            .def("calculate", 
                 py::overload_cast<const Molecule&, const MBEConfig&>(&MBECalculator::calculate),
                 py::arg("molecule"), py::arg("config"))
            
            // Analysis methods
            .def("estimate_cost", &MBECalculator::estimate_cost,
                 py::arg("molecule"), py::arg("config"))
            .def("preview_fragments", &MBECalculator::preview_fragments,
                 py::arg("molecule"), py::arg("config"))
            .def("validate_setup", &MBECalculator::validate_setup,
                 py::arg("molecule"), py::arg("config"))
            
            // Configuration
            .def("set_n_threads", &MBECalculator::set_n_threads, py::arg("n_threads"))
            .def("n_threads", &MBECalculator::n_threads)
            .def("set_caching", &MBECalculator::set_caching, py::arg("enable"))
            .def("caching_enabled", &MBECalculator::caching_enabled)
            
            // State
            .def("is_ready", &MBECalculator::is_ready)
            .def("get_status", &MBECalculator::get_status)
            .def("clear_cache", &MBECalculator::clear_cache)
            .def("get_performance_stats", &MBECalculator::get_performance_stats)
            
            .def("__repr__", [](const MBECalculator& calc) {
                return "<MBECalculator ready=" + std::string(calc.is_ready() ? "True" : "False") + ">";
            });
        
        // Convenience functions
        m.def("create_mbe_calculator", &create_mbe_calculator,
              py::arg("qm_software"), py::arg("options") = std::unordered_map<std::string, std::string>{});
        m.def("calculate_mbe", &calculate_mbe,
              py::arg("molecule"), py::arg("max_order"),
              py::arg("qm_method") = "HF", py::arg("basis_set") = "6-31G*");
    }

    void def_mbe(py::module & m) {
        // Create MBE submodule
        py::module mbe_module = m.def_submodule("mbe", "Many-Body Expansion functionality");
        
        def_mbe_config(mbe_module);
        def_mbe_results(mbe_module);
        def_mbe_fragment_generator(mbe_module);
        def_mbe_calculator(mbe_module);
        
        // Add module-level documentation
        mbe_module.doc() = R"pbdoc(
            Many-Body Expansion (MBE) module for libfrag
            
            This module provides comprehensive functionality for performing Many-Body Expansion
            calculations on molecular systems. The MBE method decomposes the total energy of
            a system into contributions from individual fragments and their interactions.
            
            Key classes:
            - MBEConfig: Configuration for MBE calculations
            - MBEResults: Storage and analysis of MBE results
            - MBEFragmentGenerator: Generation of fragment combinations
            - MBECalculator: Main calculation orchestration
        )pbdoc";
    }

} // namespace libfrag
