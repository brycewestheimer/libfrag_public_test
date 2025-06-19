#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/chrono.h>
#include <pybind11/numpy.h>
#include "xtensor-python/pyarray.hpp"

#include "libfrag/global_scf_config.hpp"
#include "libfrag/global_scf_results.hpp"
#include "libfrag/global_scf_fragment_generator.hpp"
#include "libfrag/global_scf_calculator.hpp"

namespace py = pybind11;

namespace libfrag {

    void def_global_scf_config(py::module& m) {
        py::enum_<GlobalSCFConfig::FragmentationScheme>(m, "GlobalSCFFragmentationScheme")
            .value("ATOMIC", GlobalSCFConfig::FragmentationScheme::ATOMIC)
            .value("MOLECULAR", GlobalSCFConfig::FragmentationScheme::MOLECULAR)
            .value("FUNCTIONAL_GROUP", GlobalSCFConfig::FragmentationScheme::FUNCTIONAL_GROUP)
            .value("CUSTOM", GlobalSCFConfig::FragmentationScheme::CUSTOM)
            .value("DISTANCE_BASED", GlobalSCFConfig::FragmentationScheme::DISTANCE_BASED)
            .value("FMO_LIKE", GlobalSCFConfig::FragmentationScheme::FMO_LIKE);
            
        py::enum_<GlobalSCFConfig::EmbeddingType>(m, "GlobalSCFEmbeddingType")
            .value("NONE", GlobalSCFConfig::EmbeddingType::NONE)
            .value("COULOMB", GlobalSCFConfig::EmbeddingType::COULOMB)
            .value("POLARIZABLE", GlobalSCFConfig::EmbeddingType::POLARIZABLE)
            .value("DENSITY_BASED", GlobalSCFConfig::EmbeddingType::DENSITY_BASED)
            .value("COMBINED", GlobalSCFConfig::EmbeddingType::COMBINED);
            
        py::enum_<GlobalSCFConfig::ConvergenceType>(m, "GlobalSCFConvergenceType")
            .value("ENERGY", GlobalSCFConfig::ConvergenceType::ENERGY)
            .value("DENSITY", GlobalSCFConfig::ConvergenceType::DENSITY)
            .value("ORBITAL", GlobalSCFConfig::ConvergenceType::ORBITAL)
            .value("COMBINED", GlobalSCFConfig::ConvergenceType::COMBINED);

        py::class_<GlobalSCFConfig>(m, "GlobalSCFConfig", R"pbdoc(
            Configuration class for Global SCF calculations.
            
            This class stores all parameters needed to configure Global SCF calculations,
            including fragmentation schemes, embedding potential types, convergence criteria,
            and quantum chemistry method specifications for fragment-based SCF methods.
            
            Examples
            --------
            Create a basic Global SCF configuration:
            
            >>> config = GlobalSCFConfig()
            >>> config.set_fragmentation_scheme(GlobalSCFFragmentationScheme.MOLECULAR)
            >>> config.set_embedding_type(GlobalSCFEmbeddingType.COULOMB)
            >>> config.set_qm_method("HF")
            >>> config.set_basis_set("6-31G*")
            
            Create FMO-style configuration:
            
            >>> fmo_config = GlobalSCFConfig.create_fmo_config("B3LYP", "def2-TZVP")
            >>> fmo_config.set_max_scf_iterations(100)
            >>> fmo_config.set_energy_threshold(1e-6)
            
            Configure for polarizable embedding:
            
            >>> pol_config = GlobalSCFConfig.create_polarizable_config()
            >>> pol_config.set_polarization_threshold(1e-6)
            >>> pol_config.set_diis_enabled(True)
        )pbdoc")
            .def(py::init<>())
            .def(py::init<GlobalSCFConfig::FragmentationScheme, GlobalSCFConfig::EmbeddingType>(),
                 py::arg("scheme"), py::arg("embedding") = GlobalSCFConfig::EmbeddingType::COULOMB)
            
            // Configuration setters
            .def("set_fragmentation_scheme", &GlobalSCFConfig::set_fragmentation_scheme, py::arg("scheme"))
            .def("set_embedding_type", &GlobalSCFConfig::set_embedding_type, py::arg("embedding"))
            .def("set_convergence_type", &GlobalSCFConfig::set_convergence_type, py::arg("convergence"))
            .def("set_max_scf_iterations", &GlobalSCFConfig::set_max_scf_iterations, py::arg("max_iter"))
            .def("set_energy_threshold", &GlobalSCFConfig::set_energy_threshold, py::arg("threshold"))
            .def("set_density_threshold", &GlobalSCFConfig::set_density_threshold, py::arg("threshold"))
            .def("set_orbital_threshold", &GlobalSCFConfig::set_orbital_threshold, py::arg("threshold"))
            .def("set_qm_method", &GlobalSCFConfig::set_qm_method, py::arg("method"))
            .def("set_basis_set", &GlobalSCFConfig::set_basis_set, py::arg("basis"))
            .def("set_diis_enabled", &GlobalSCFConfig::set_diis_enabled, py::arg("enable"))
            .def("set_diis_size", &GlobalSCFConfig::set_diis_size, py::arg("size"))
            .def("set_buffer_distance", &GlobalSCFConfig::set_buffer_distance, py::arg("distance"))
            .def("set_custom_fragments", &GlobalSCFConfig::set_custom_fragments, py::arg("fragments"))
            .def("set_qm_options", &GlobalSCFConfig::set_qm_options, py::arg("options"))
            .def("set_polarization_enabled", &GlobalSCFConfig::set_polarization_enabled, py::arg("enable"))
            .def("set_polarization_threshold", &GlobalSCFConfig::set_polarization_threshold, py::arg("threshold"))
            
            // Configuration getters
            .def("fragmentation_scheme", &GlobalSCFConfig::fragmentation_scheme)
            .def("embedding_type", &GlobalSCFConfig::embedding_type)
            .def("convergence_type", &GlobalSCFConfig::convergence_type)
            .def("max_scf_iterations", &GlobalSCFConfig::max_scf_iterations)
            .def("energy_threshold", &GlobalSCFConfig::energy_threshold)
            .def("density_threshold", &GlobalSCFConfig::density_threshold)
            .def("orbital_threshold", &GlobalSCFConfig::orbital_threshold)
            .def("qm_method", &GlobalSCFConfig::qm_method)
            .def("basis_set", &GlobalSCFConfig::basis_set)
            .def("diis_enabled", &GlobalSCFConfig::diis_enabled)
            .def("diis_size", &GlobalSCFConfig::diis_size)
            .def("buffer_distance", &GlobalSCFConfig::buffer_distance)
            .def("polarization_enabled", &GlobalSCFConfig::polarization_enabled)
            .def("polarization_threshold", &GlobalSCFConfig::polarization_threshold)
            .def("custom_fragments", &GlobalSCFConfig::custom_fragments)
            .def("qm_options", &GlobalSCFConfig::qm_options)
            
            // Validation and utilities
            .def("is_valid", &GlobalSCFConfig::is_valid)
            .def("to_string", &GlobalSCFConfig::to_string)
            
            // Factory methods
            .def_static("create_fmo_config", &GlobalSCFConfig::create_fmo_config,
                       py::arg("method") = "HF", py::arg("basis") = "6-31G*")
            .def_static("create_coulomb_config", &GlobalSCFConfig::create_coulomb_config,
                       py::arg("method") = "HF", py::arg("basis") = "6-31G*")
            .def_static("create_polarizable_config", &GlobalSCFConfig::create_polarizable_config,
                       py::arg("method") = "HF", py::arg("basis") = "6-31G*")
            
            // String conversion utilities
            .def_static("fragmentation_scheme_to_string", &GlobalSCFConfig::fragmentation_scheme_to_string)
            .def_static("embedding_type_to_string", &GlobalSCFConfig::embedding_type_to_string)
            .def_static("convergence_type_to_string", &GlobalSCFConfig::convergence_type_to_string);
    }

    void def_global_scf_results(py::module& m) {
        py::class_<FragmentWaveFunction>(m, "FragmentWaveFunction", R"pbdoc(
            Fragment wave function data container.
            
            Stores the quantum chemistry calculation results and wave function data
            for a single fragment in the Global SCF calculation.
        )pbdoc")
            .def(py::init<>())
            .def(py::init<std::size_t, const std::vector<std::size_t>&>(),
                 py::arg("fragment_id"), py::arg("atom_indices"))
            
            // Data members
            .def_readwrite("fragment_id", &FragmentWaveFunction::fragment_id)
            .def_readwrite("atom_indices", &FragmentWaveFunction::atom_indices)
            .def_readwrite("fragment_name", &FragmentWaveFunction::fragment_name)
            .def_readwrite("total_energy", &FragmentWaveFunction::total_energy)
            .def_readwrite("electronic_energy", &FragmentWaveFunction::electronic_energy)
            .def_readwrite("nuclear_repulsion", &FragmentWaveFunction::nuclear_repulsion)
            .def_readwrite("embedding_energy", &FragmentWaveFunction::embedding_energy)
            .def_readwrite("polarization_energy", &FragmentWaveFunction::polarization_energy)
            .def_readwrite("qm_method", &FragmentWaveFunction::qm_method)
            .def_readwrite("basis_set", &FragmentWaveFunction::basis_set)
            .def_readwrite("n_electrons", &FragmentWaveFunction::n_electrons)
            .def_readwrite("multiplicity", &FragmentWaveFunction::multiplicity)
            .def_readwrite("charge", &FragmentWaveFunction::charge)
            .def_readwrite("converged", &FragmentWaveFunction::converged)
            .def_readwrite("scf_iterations", &FragmentWaveFunction::scf_iterations)
            .def_readwrite("computation_time", &FragmentWaveFunction::computation_time)
            .def_readwrite("status", &FragmentWaveFunction::status)
            .def_readwrite("properties", &FragmentWaveFunction::properties)
            
            // Methods
            .def("to_string", &FragmentWaveFunction::to_string)
            .def("is_valid", &FragmentWaveFunction::is_valid)
            .def("get_property", &FragmentWaveFunction::get_property, 
                 py::arg("name"), py::arg("default_value") = 0.0)
            .def("set_property", &FragmentWaveFunction::set_property, 
                 py::arg("name"), py::arg("value"))
            .def("electronic_charge", &FragmentWaveFunction::electronic_charge)
            .def("nuclear_charge", &FragmentWaveFunction::nuclear_charge)
            .def("electric_dipole_moment", &FragmentWaveFunction::electric_dipole_moment)
            .def("total_electron_density", &FragmentWaveFunction::total_electron_density);

        py::class_<GlobalSCFIteration>(m, "GlobalSCFIteration", R"pbdoc(
            Global SCF iteration data.
            
            Stores information about a single global SCF iteration,
            including energies, convergence criteria, and timing.
        )pbdoc")
            .def(py::init<>())
            .def(py::init<int>(), py::arg("iteration_number"))
            
            // Data members
            .def_readwrite("iteration_number", &GlobalSCFIteration::iteration_number)
            .def_readwrite("total_energy", &GlobalSCFIteration::total_energy)
            .def_readwrite("electronic_energy", &GlobalSCFIteration::electronic_energy)
            .def_readwrite("nuclear_repulsion", &GlobalSCFIteration::nuclear_repulsion)
            .def_readwrite("embedding_energy", &GlobalSCFIteration::embedding_energy)
            .def_readwrite("polarization_energy", &GlobalSCFIteration::polarization_energy)
            .def_readwrite("interaction_energy", &GlobalSCFIteration::interaction_energy)
            .def_readwrite("energy_change", &GlobalSCFIteration::energy_change)
            .def_readwrite("max_density_change", &GlobalSCFIteration::max_density_change)
            .def_readwrite("rms_density_change", &GlobalSCFIteration::rms_density_change)
            .def_readwrite("max_orbital_change", &GlobalSCFIteration::max_orbital_change)
            .def_readwrite("polarization_change", &GlobalSCFIteration::polarization_change)
            .def_readwrite("iteration_time", &GlobalSCFIteration::iteration_time)
            .def_readwrite("fragment_time", &GlobalSCFIteration::fragment_time)
            .def_readwrite("embedding_time", &GlobalSCFIteration::embedding_time)
            .def_readwrite("energy_converged", &GlobalSCFIteration::energy_converged)
            .def_readwrite("density_converged", &GlobalSCFIteration::density_converged)
            .def_readwrite("orbital_converged", &GlobalSCFIteration::orbital_converged)
            .def_readwrite("polarization_converged", &GlobalSCFIteration::polarization_converged)
            .def_readwrite("overall_converged", &GlobalSCFIteration::overall_converged)
            .def_readwrite("diis_used", &GlobalSCFIteration::diis_used)
            .def_readwrite("diis_dimension", &GlobalSCFIteration::diis_dimension)
            .def_readwrite("diis_error", &GlobalSCFIteration::diis_error)
            .def_readwrite("status_message", &GlobalSCFIteration::status_message)
            .def_readwrite("diagnostics", &GlobalSCFIteration::diagnostics)
            
            // Methods
            .def("to_string", &GlobalSCFIteration::to_string)
            .def("is_converged", &GlobalSCFIteration::is_converged, py::arg("config"))
            .def("set_diagnostic", &GlobalSCFIteration::set_diagnostic, 
                 py::arg("name"), py::arg("value"))
            .def("get_diagnostic", &GlobalSCFIteration::get_diagnostic, 
                 py::arg("name"), py::arg("default_value") = 0.0);

        py::class_<GlobalSCFResults>(m, "GlobalSCFResults", R"pbdoc(
            Results container for Global SCF calculations.
            
            This class stores the complete results of a Global SCF calculation,
            including fragment wave functions, iteration history, convergence analysis,
            and performance statistics.
            
            Examples
            --------
            Access calculation results:
            
            >>> results = calculator.calculate(molecule)
            >>> print(f"Final energy: {results.final_energy()} Eh")
            >>> print(f"Converged: {results.is_converged()}")
            >>> print(f"Iterations: {results.n_iterations()}")
            
            Analyze convergence:
            
            >>> conv_analysis = results.convergence_analysis()
            >>> energy_history = results.convergence_history("energy_change")
            >>> density_history = results.convergence_history("density_change")
            
            Export results:
            
            >>> results.export_json("global_scf_results.json")
            >>> results.export_convergence_csv("convergence.csv")
            >>> print(results.generate_summary(detailed=True))
        )pbdoc")
            .def(py::init<>())
            .def(py::init<const GlobalSCFConfig&>(), py::arg("config"))
            
            // Fragment wave function management
            .def("add_fragment_wavefunction", &GlobalSCFResults::add_fragment_wavefunction, py::arg("wf"))
            .def("get_fragment_wavefunction", 
                 py::overload_cast<std::size_t>(&GlobalSCFResults::get_fragment_wavefunction, py::const_),
                 py::arg("fragment_id"), py::return_value_policy::reference_internal)
            .def("fragment_wavefunctions", &GlobalSCFResults::fragment_wavefunctions,
                 py::return_value_policy::reference_internal)
            
            // Iteration history management
            .def("add_scf_iteration", &GlobalSCFResults::add_scf_iteration, py::arg("iteration"))
            .def("get_scf_iteration", &GlobalSCFResults::get_scf_iteration, py::arg("iteration_number"),
                 py::return_value_policy::reference_internal)
            .def("scf_iterations", &GlobalSCFResults::scf_iterations,
                 py::return_value_policy::reference_internal)
            
            // Energy analysis
            .def("final_energy", &GlobalSCFResults::final_energy)
            .def("energy_breakdown", &GlobalSCFResults::energy_breakdown)
            .def("fragment_energies", &GlobalSCFResults::fragment_energies)
            .def("interaction_energy", &GlobalSCFResults::interaction_energy)
            
            // Convergence analysis
            .def("is_converged", &GlobalSCFResults::is_converged)
            .def("n_iterations", &GlobalSCFResults::n_iterations)
            .def("convergence_analysis", &GlobalSCFResults::convergence_analysis)
            .def("convergence_history", &GlobalSCFResults::convergence_history, py::arg("criterion"))
            
            // Performance analysis
            .def("total_time", &GlobalSCFResults::total_time)
            .def("performance_breakdown", &GlobalSCFResults::performance_breakdown)
            .def("performance_statistics", &GlobalSCFResults::performance_statistics)
            
            // Export and analysis
            .def("export_json", &GlobalSCFResults::export_json, 
                 py::arg("filename"), py::arg("include_wavefunctions") = false)
            .def("export_convergence_csv", &GlobalSCFResults::export_convergence_csv, py::arg("filename"))
            .def("export_energy_csv", &GlobalSCFResults::export_energy_csv, py::arg("filename"))
            .def("generate_summary", &GlobalSCFResults::generate_summary, py::arg("detailed") = false)
            
            // Validation and diagnostics
            .def("validate", &GlobalSCFResults::validate)
            .def("calculation_statistics", &GlobalSCFResults::calculation_statistics)
            
            // Configuration and utility
            .def("config", &GlobalSCFResults::config, py::return_value_policy::reference_internal)
            .def("set_config", &GlobalSCFResults::set_config, py::arg("config"))
            .def("n_fragments", &GlobalSCFResults::n_fragments)
            .def("empty", &GlobalSCFResults::empty)
            .def("clear", &GlobalSCFResults::clear)
            .def("to_string", &GlobalSCFResults::to_string);
    }

    void def_global_scf_fragment_generator(py::module& m) {
        py::class_<GlobalSCFFragment>(m, "GlobalSCFFragment", R"pbdoc(
            Fragment definition for Global SCF calculations.
            
            Extended fragment information including buffer regions and
            embedding environment for Global SCF calculations.
        )pbdoc")
            .def(py::init<>())
            .def(py::init<std::size_t, std::shared_ptr<Fragment>>(),
                 py::arg("fragment_id"), py::arg("core_fragment"))
            
            // Data members
            .def_readwrite("fragment_id", &GlobalSCFFragment::fragment_id)
            .def_readwrite("fragment_name", &GlobalSCFFragment::fragment_name)
            .def_readwrite("core_atoms", &GlobalSCFFragment::core_atoms)
            .def_readwrite("buffer_atoms", &GlobalSCFFragment::buffer_atoms)
            .def_readwrite("environment_atoms", &GlobalSCFFragment::environment_atoms)
            .def_readwrite("center_of_mass", &GlobalSCFFragment::center_of_mass)
            .def_readwrite("radius", &GlobalSCFFragment::radius)
            .def_readwrite("bounding_box", &GlobalSCFFragment::bounding_box)
            .def_readwrite("neighboring_fragments", &GlobalSCFFragment::neighboring_fragments)
            .def_readwrite("fragment_distances", &GlobalSCFFragment::fragment_distances)
            .def_readwrite("buffer_distance", &GlobalSCFFragment::buffer_distance)
            .def_readwrite("capping_atoms", &GlobalSCFFragment::capping_atoms)
            .def_readwrite("total_charge", &GlobalSCFFragment::total_charge)
            .def_readwrite("multiplicity", &GlobalSCFFragment::multiplicity)
            
            // Methods
            .def("to_string", &GlobalSCFFragment::to_string)
            .def("is_valid", &GlobalSCFFragment::is_valid)
            .def("total_atoms", &GlobalSCFFragment::total_atoms)
            .def("contains_atom", &GlobalSCFFragment::contains_atom, py::arg("atom_index"))
            .def("is_neighbor", &GlobalSCFFragment::is_neighbor, py::arg("other_fragment_id"))
            .def("distance_to_fragment", &GlobalSCFFragment::distance_to_fragment, py::arg("other_fragment_id"));

        py::class_<GlobalSCFFragmentGenerator>(m, "GlobalSCFFragmentGenerator", R"pbdoc(
            Generates fragment definitions for Global SCF calculations.
            
            This class is responsible for creating fragment definitions suitable for
            Global SCF calculations, including core fragments, buffer regions, and
            embedding environments.
            
            Examples
            --------
            Generate fragments for a molecule:
            
            >>> generator = GlobalSCFFragmentGenerator()
            >>> config = GlobalSCFConfig.create_fmo_config()
            >>> fragments = generator.generate_fragments(molecule, config)
            
            Analyze fragmentation:
            
            >>> stats = generator.get_fragmentation_statistics(fragments)
            >>> cost = generator.estimate_computational_cost(fragments, config)
            >>> report = generator.generate_report(fragments, molecule)
        )pbdoc")
            .def(py::init<>())
            .def(py::init<const GlobalSCFConfig&>(), py::arg("config"))
            
            // Configuration
            .def("set_config", &GlobalSCFFragmentGenerator::set_config, py::arg("config"))
            .def("config", &GlobalSCFFragmentGenerator::config, py::return_value_policy::reference_internal)
            
            // Main fragmentation methods
            .def("generate_fragments", 
                 py::overload_cast<const Molecule&>(&GlobalSCFFragmentGenerator::generate_fragments),
                 py::arg("molecule"))
            .def("generate_fragments", 
                 py::overload_cast<const Molecule&, const GlobalSCFConfig&>(&GlobalSCFFragmentGenerator::generate_fragments),
                 py::arg("molecule"), py::arg("config"))
            
            // Specific fragmentation schemes
            .def("fragment_by_atoms", &GlobalSCFFragmentGenerator::fragment_by_atoms, py::arg("molecule"))
            .def("fragment_by_molecules", &GlobalSCFFragmentGenerator::fragment_by_molecules, py::arg("molecule"))
            .def("fragment_by_functional_groups", &GlobalSCFFragmentGenerator::fragment_by_functional_groups, py::arg("molecule"))
            .def("fragment_by_distance", &GlobalSCFFragmentGenerator::fragment_by_distance,
                 py::arg("molecule"), py::arg("cutoff"))
            .def("fragment_by_custom", &GlobalSCFFragmentGenerator::fragment_by_custom,
                 py::arg("molecule"), py::arg("fragment_definitions"))
            .def("fragment_fmo_style", &GlobalSCFFragmentGenerator::fragment_fmo_style, py::arg("molecule"))
            
            // Fragment analysis
            .def("determine_neighbors", &GlobalSCFFragmentGenerator::determine_neighbors,
                 py::arg("fragments"), py::arg("molecule"))
            .def("calculate_distances", &GlobalSCFFragmentGenerator::calculate_distances, py::arg("fragments"))
            .def("analyze_overlap", &GlobalSCFFragmentGenerator::analyze_overlap, py::arg("fragments"))
            
            // Validation and optimization
            .def("validate_fragments", &GlobalSCFFragmentGenerator::validate_fragments,
                 py::arg("fragments"), py::arg("molecule"))
            .def("optimize_boundaries", &GlobalSCFFragmentGenerator::optimize_boundaries,
                 py::arg("fragments"), py::arg("molecule"))
            .def("balance_fragment_sizes", &GlobalSCFFragmentGenerator::balance_fragment_sizes,
                 py::arg("fragments"), py::arg("molecule"), py::arg("target_size"))
            
            // Statistics and analysis
            .def("get_fragmentation_statistics", &GlobalSCFFragmentGenerator::get_fragmentation_statistics, py::arg("fragments"))
            .def("estimate_computational_cost", &GlobalSCFFragmentGenerator::estimate_computational_cost,
                 py::arg("fragments"), py::arg("config"))
            .def("generate_report", &GlobalSCFFragmentGenerator::generate_report,
                 py::arg("fragments"), py::arg("molecule"))
            
            // Utility methods
            .def("create_fragment_from_atoms", &GlobalSCFFragmentGenerator::create_fragment_from_atoms,
                 py::arg("molecule"), py::arg("atom_indices"), py::arg("fragment_id"))
            .def("find_atoms_within_distance", &GlobalSCFFragmentGenerator::find_atoms_within_distance,
                 py::arg("molecule"), py::arg("fragment"), py::arg("distance"))
            .def("get_fragment_containing_atom", &GlobalSCFFragmentGenerator::get_fragment_containing_atom,
                 py::arg("fragments"), py::arg("atom_index"));
    }

    void def_global_scf_calculator(py::module& m) {
        // QM Interface (abstract base class)
        py::class_<GlobalSCFQMInterface>(m, "GlobalSCFQMInterface", R"pbdoc(
            Interface for quantum chemistry calculations in Global SCF.
            
            Abstract base class defining interface for quantum chemistry calculations
            needed by the Global SCF calculator. This extends the basic QM interface
            to support embedding potentials and polarization effects.
        )pbdoc")
            .def("calculate_embedded_scf", &GlobalSCFQMInterface::calculate_embedded_scf,
                 py::arg("fragment"), py::arg("method"), py::arg("basis_set"),
                 py::arg("embedding_potential"), py::arg("point_charges"), py::arg("charge_positions"),
                 py::arg("charge") = 0, py::arg("multiplicity") = 1,
                 py::arg("options") = std::unordered_map<std::string, std::string>{})
            .def("calculate_in_vacuo", &GlobalSCFQMInterface::calculate_in_vacuo,
                 py::arg("fragment"), py::arg("method"), py::arg("basis_set"),
                 py::arg("charge") = 0, py::arg("multiplicity") = 1,
                 py::arg("options") = std::unordered_map<std::string, std::string>{})
            .def("calculate_polarization_response", &GlobalSCFQMInterface::calculate_polarization_response,
                 py::arg("fragment"), py::arg("external_field"), py::arg("method"), py::arg("basis_set"),
                 py::arg("options") = std::unordered_map<std::string, std::string>{})
            .def("software_info", &GlobalSCFQMInterface::software_info)
            .def("is_available", &GlobalSCFQMInterface::is_available)
            .def("supported_methods", &GlobalSCFQMInterface::supported_methods)
            .def("supported_basis_sets", &GlobalSCFQMInterface::supported_basis_sets);

        py::class_<DIISAccelerator>(m, "DIISAccelerator", R"pbdoc(
            DIIS (Direct Inversion of Iterative Subspace) accelerator.
            
            Implements DIIS acceleration for Global SCF convergence,
            working with density matrices and/or Fock matrices.
        )pbdoc")
            .def(py::init<int>(), py::arg("max_vectors") = 8)
            .def("add_iteration", &DIISAccelerator::add_iteration,
                 py::arg("density"), py::arg("fock"), py::arg("error"))
            .def("get_extrapolated_density", &DIISAccelerator::get_extrapolated_density)
            .def("get_extrapolated_fock", &DIISAccelerator::get_extrapolated_fock)
            .def("is_ready", &DIISAccelerator::is_ready)
            .def("current_error", &DIISAccelerator::current_error)
            .def("reset", &DIISAccelerator::reset)
            .def("subspace_dimension", &DIISAccelerator::subspace_dimension);

        py::class_<GlobalSCFCalculator>(m, "GlobalSCFCalculator", R"pbdoc(
            Main Global SCF calculator.
            
            This class orchestrates the complete Global SCF calculation workflow:
            1. Fragment generation and setup
            2. Initial fragment calculations (in vacuo or with simple embedding)
            3. Embedding potential generation from fragment wave functions
            4. Iterative fragment recalculation with updated embedding
            5. Convergence checking and DIIS acceleration
            6. Result compilation and analysis
            
            Examples
            --------
            Basic Global SCF calculation:
            
            >>> calculator = GlobalSCFCalculator()
            >>> config = GlobalSCFConfig.create_fmo_config("HF", "6-31G*")
            >>> results = calculator.calculate(molecule, config)
            
            With progress monitoring:
            
            >>> def progress_callback(iteration, progress, status):
            ...     print(f"Iteration {iteration}: {progress:.1f}% - {status}")
            >>> calculator.set_progress_callback(progress_callback)
            >>> results = calculator.calculate(molecule)
            
            Async calculation:
            
            >>> future_results = calculator.calculate_async(molecule)
            >>> # Do other work...
            >>> results = future_results.get()
        )pbdoc")
            .def(py::init<>())
            .def(py::init<std::unique_ptr<GlobalSCFQMInterface>>(), py::arg("qm_interface"))
            .def(py::init<std::unique_ptr<GlobalSCFQMInterface>, const GlobalSCFConfig&>(),
                 py::arg("qm_interface"), py::arg("config"))
            
            // Configuration
            .def("set_config", &GlobalSCFCalculator::set_config, py::arg("config"))
            .def("config", &GlobalSCFCalculator::config, py::return_value_policy::reference_internal)
            .def("set_qm_interface", &GlobalSCFCalculator::set_qm_interface, py::arg("qm_interface"))
            
            // Callback setup
            .def("set_progress_callback", &GlobalSCFCalculator::set_progress_callback, py::arg("callback"))
            .def("set_log_callback", &GlobalSCFCalculator::set_log_callback, py::arg("callback"))
            .def("set_convergence_callback", &GlobalSCFCalculator::set_convergence_callback, py::arg("callback"))
            
            // Main calculation methods
            .def("calculate", 
                 py::overload_cast<const Molecule&>(&GlobalSCFCalculator::calculate),
                 py::arg("molecule"))
            .def("calculate", 
                 py::overload_cast<const Molecule&, const GlobalSCFConfig&>(&GlobalSCFCalculator::calculate),
                 py::arg("molecule"), py::arg("config"))
            .def("calculate_async", &GlobalSCFCalculator::calculate_async, py::arg("molecule"))
            
            // Incremental calculation methods
            .def("start_calculation", &GlobalSCFCalculator::start_calculation, py::arg("molecule"))
            .def("get_partial_results", &GlobalSCFCalculator::get_partial_results, py::arg("calculation_id"))
            .def("is_calculation_complete", &GlobalSCFCalculator::is_calculation_complete, py::arg("calculation_id"))
            .def("cancel_calculation", &GlobalSCFCalculator::cancel_calculation, py::arg("calculation_id"))
            
            // Analysis and prediction methods
            .def("estimate_cost", &GlobalSCFCalculator::estimate_cost,
                 py::arg("molecule"), py::arg("config"))
            .def("preview_fragments", &GlobalSCFCalculator::preview_fragments,
                 py::arg("molecule"), py::arg("config"))
            .def("validate_setup", &GlobalSCFCalculator::validate_setup,
                 py::arg("molecule"), py::arg("config"))
            
            // Parallelization control
            .def("set_n_threads", &GlobalSCFCalculator::set_n_threads, py::arg("n_threads"))
            .def("n_threads", &GlobalSCFCalculator::n_threads)
            .def("set_caching", &GlobalSCFCalculator::set_caching, py::arg("enable"))
            .def("caching_enabled", &GlobalSCFCalculator::caching_enabled)
            
            // State and diagnostics
            .def("is_ready", &GlobalSCFCalculator::is_ready)
            .def("get_status", &GlobalSCFCalculator::get_status)
            .def("get_performance_stats", &GlobalSCFCalculator::get_performance_stats);

        // Factory and utility functions
        m.def("create_global_scf_calculator", &create_global_scf_calculator,
              py::arg("qm_software"), py::arg("options") = std::unordered_map<std::string, std::string>{});
        m.def("calculate_global_scf", &calculate_global_scf,
              py::arg("molecule"), 
              py::arg("embedding_type") = GlobalSCFConfig::EmbeddingType::COULOMB,
              py::arg("qm_method") = "HF", py::arg("basis_set") = "6-31G*");
    }

    void def_global_scf(py::module& m) {
        def_global_scf_config(m);
        def_global_scf_results(m);
        def_global_scf_fragment_generator(m);
        def_global_scf_calculator(m);
    }

}
