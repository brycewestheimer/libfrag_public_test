#pragma once
#ifndef LIBFRAG_GLOBAL_SCF_CALCULATOR_HPP
#define LIBFRAG_GLOBAL_SCF_CALCULATOR_HPP

#include <memory>
#include <vector>
#include <unordered_map>
#include <string>
#include <future>
#include <functional>
#include <chrono>
#include "molecule.hpp"
#include "fragment.hpp"
#include "global_scf_config.hpp"
#include "global_scf_results.hpp"
#include "global_scf_fragment_generator.hpp"
#include "xtensor/xarray.hpp"

namespace libfrag {

    /**
     * @brief Interface for quantum chemistry calculations in Global SCF
     * 
     * Abstract base class defining interface for quantum chemistry calculations
     * needed by the Global SCF calculator. This extends the basic QM interface
     * to support embedding potentials and polarization effects.
     */
    class GlobalSCFQMInterface {
    public:
        virtual ~GlobalSCFQMInterface() = default;

        /**
         * @brief Perform SCF calculation with embedding potential
         * @param fragment Fragment to calculate
         * @param method QM method string
         * @param basis_set Basis set string
         * @param embedding_potential External embedding potential
         * @param point_charges Point charges from environment
         * @param charge_positions Positions of point charges
         * @param charge Total charge
         * @param multiplicity Spin multiplicity
         * @param options Additional QM options
         * @return Fragment wave function result
         */
        virtual FragmentWaveFunction calculate_embedded_scf(
            const Fragment& fragment,
            const std::string& method,
            const std::string& basis_set,
            const xt::xarray<double>& embedding_potential,
            const std::vector<double>& point_charges,
            const std::vector<std::array<double, 3>>& charge_positions,
            int charge = 0,
            int multiplicity = 1,
            const std::unordered_map<std::string, std::string>& options = {}) = 0;

        /**
         * @brief Calculate in vacuo (no embedding)
         * @param fragment Fragment to calculate
         * @param method QM method string
         * @param basis_set Basis set string
         * @param charge Total charge
         * @param multiplicity Spin multiplicity
         * @param options Additional QM options
         * @return Fragment wave function result
         */
        virtual FragmentWaveFunction calculate_in_vacuo(
            const Fragment& fragment,
            const std::string& method,
            const std::string& basis_set,
            int charge = 0,
            int multiplicity = 1,
            const std::unordered_map<std::string, std::string>& options = {}) = 0;

        /**
         * @brief Calculate polarization response
         * @param fragment Fragment to calculate
         * @param external_field External electric field
         * @param method QM method string
         * @param basis_set Basis set string
         * @param options Additional QM options
         * @return Polarization response data
         */
        virtual xt::xarray<double> calculate_polarization_response(
            const Fragment& fragment,
            const xt::xarray<double>& external_field,
            const std::string& method,
            const std::string& basis_set,
            const std::unordered_map<std::string, std::string>& options = {}) = 0;

        /**
         * @brief Get software information
         * @return Software name and version
         */
        virtual std::string software_info() const = 0;

        /**
         * @brief Check if software is available
         * @return True if software is available
         */
        virtual bool is_available() const = 0;

        /**
         * @brief Get supported methods
         * @return List of supported QM methods
         */
        virtual std::vector<std::string> supported_methods() const = 0;

        /**
         * @brief Get supported basis sets
         * @return List of supported basis sets
         */
        virtual std::vector<std::string> supported_basis_sets() const = 0;
    };

    /**
     * @brief DIIS (Direct Inversion of Iterative Subspace) accelerator
     * 
     * Implements DIIS acceleration for Global SCF convergence,
     * working with density matrices and/or Fock matrices.
     */
    class DIISAccelerator {
    public:
        /**
         * @brief Construct DIIS accelerator
         * @param max_vectors Maximum number of DIIS vectors
         */
        explicit DIISAccelerator(int max_vectors = 8);

        /**
         * @brief Add iteration data to DIIS
         * @param density Current density matrix
         * @param fock Current Fock matrix
         * @param error Error vector for this iteration
         */
        void add_iteration(const xt::xarray<double>& density,
                          const xt::xarray<double>& fock,
                          const xt::xarray<double>& error);

        /**
         * @brief Get DIIS-extrapolated density matrix
         * @return Extrapolated density matrix
         */
        xt::xarray<double> get_extrapolated_density();

        /**
         * @brief Get DIIS-extrapolated Fock matrix
         * @return Extrapolated Fock matrix
         */
        xt::xarray<double> get_extrapolated_fock();

        /**
         * @brief Check if DIIS is ready (has enough vectors)
         * @return True if DIIS can be used
         */
        bool is_ready() const { return stored_vectors_ >= 2; }

        /**
         * @brief Get current DIIS error
         * @return DIIS error norm
         */
        double current_error() const { return current_error_; }

        /**
         * @brief Reset DIIS (clear all stored vectors)
         */
        void reset();

        /**
         * @brief Get DIIS subspace dimension
         * @return Number of stored vectors
         */
        int subspace_dimension() const { return stored_vectors_; }

    private:
        int max_vectors_;
        int stored_vectors_ = 0;
        int current_index_ = 0;
        double current_error_ = 0.0;

        std::vector<xt::xarray<double>> density_vectors_;
        std::vector<xt::xarray<double>> fock_vectors_;
        std::vector<xt::xarray<double>> error_vectors_;

        // DIIS coefficient calculation
        std::vector<double> calculate_diis_coefficients();
        double calculate_error_norm(const xt::xarray<double>& error);
    };

    /**
     * @brief Main Global SCF calculator
     * 
     * This class orchestrates the complete Global SCF calculation workflow:
     * 1. Fragment generation and setup
     * 2. Initial fragment calculations (in vacuo or with simple embedding)
     * 3. Embedding potential generation from fragment wave functions
     * 4. Iterative fragment recalculation with updated embedding
     * 5. Convergence checking and DIIS acceleration
     * 6. Result compilation and analysis
     */
    class GlobalSCFCalculator {
    public:
        // Type aliases
        using QMInterfacePtr = std::unique_ptr<GlobalSCFQMInterface>;
        using ProgressCallback = std::function<void(int iteration, double progress, const std::string& status)>;
        using LogCallback = std::function<void(const std::string& message)>;
        using ConvergenceCallback = std::function<void(const GlobalSCFIteration& iteration)>;

        // Constructors
        GlobalSCFCalculator() = default;
        
        /**
         * @brief Construct with QM interface
         * @param qm_interface Quantum chemistry interface implementation
         */
        explicit GlobalSCFCalculator(QMInterfacePtr qm_interface);

        /**
         * @brief Construct with QM interface and configuration
         * @param qm_interface Quantum chemistry interface implementation
         * @param config Global SCF configuration
         */
        GlobalSCFCalculator(QMInterfacePtr qm_interface, const GlobalSCFConfig& config);

        // Configuration
        /**
         * @brief Set Global SCF configuration
         * @param config Global SCF configuration
         */
        void set_config(const GlobalSCFConfig& config) { config_ = config; }

        /**
         * @brief Get current configuration
         * @return Current configuration
         */
        const GlobalSCFConfig& config() const { return config_; }

        /**
         * @brief Set quantum chemistry interface
         * @param qm_interface QM interface implementation
         */
        void set_qm_interface(QMInterfacePtr qm_interface);

        // Callback setup
        /**
         * @brief Set progress callback
         * @param callback Progress callback function
         */
        void set_progress_callback(ProgressCallback callback) { progress_callback_ = callback; }

        /**
         * @brief Set logging callback
         * @param callback Log callback function
         */
        void set_log_callback(LogCallback callback) { log_callback_ = callback; }

        /**
         * @brief Set convergence callback
         * @param callback Convergence callback function
         */
        void set_convergence_callback(ConvergenceCallback callback) { convergence_callback_ = callback; }

        // Main calculation methods
        /**
         * @brief Perform Global SCF calculation
         * @param molecule Input molecule
         * @return Global SCF calculation results
         */
        GlobalSCFResults calculate(const Molecule& molecule);

        /**
         * @brief Perform Global SCF calculation with custom configuration
         * @param molecule Input molecule
         * @param config Global SCF configuration
         * @return Global SCF calculation results
         */
        GlobalSCFResults calculate(const Molecule& molecule, const GlobalSCFConfig& config);

        /**
         * @brief Perform async Global SCF calculation
         * @param molecule Input molecule
         * @return Future for Global SCF calculation results
         */
        std::future<GlobalSCFResults> calculate_async(const Molecule& molecule);

        // Incremental calculation methods
        /**
         * @brief Start incremental Global SCF calculation
         * @param molecule Input molecule
         * @return Calculation handle for incremental updates
         */
        std::string start_calculation(const Molecule& molecule);

        /**
         * @brief Get partial results from ongoing calculation
         * @param calculation_id Calculation identifier
         * @return Partial results
         */
        GlobalSCFResults get_partial_results(const std::string& calculation_id);

        /**
         * @brief Check if calculation is complete
         * @param calculation_id Calculation identifier
         * @return True if calculation is finished
         */
        bool is_calculation_complete(const std::string& calculation_id);

        /**
         * @brief Cancel ongoing calculation
         * @param calculation_id Calculation identifier
         */
        void cancel_calculation(const std::string& calculation_id);

        // Analysis and prediction methods
        /**
         * @brief Estimate computational cost
         * @param molecule Input molecule
         * @param config Global SCF configuration
         * @return Estimated cost metrics
         */
        std::unordered_map<std::string, double> estimate_cost(
            const Molecule& molecule, const GlobalSCFConfig& config);

        /**
         * @brief Preview fragment definitions
         * @param molecule Input molecule
         * @param config Global SCF configuration
         * @return Fragment preview information
         */
        std::unordered_map<std::string, double> preview_fragments(
            const Molecule& molecule, const GlobalSCFConfig& config);

        /**
         * @brief Validate calculation setup
         * @param molecule Input molecule
         * @param config Global SCF configuration
         * @return Validation report
         */
        std::vector<std::string> validate_setup(
            const Molecule& molecule, const GlobalSCFConfig& config);

        // Parallelization control
        /**
         * @brief Set number of parallel threads for fragment calculations
         * @param n_threads Number of threads (0 = auto-detect)
         */
        void set_n_threads(int n_threads) { n_threads_ = n_threads; }

        /**
         * @brief Get number of parallel threads
         * @return Number of threads being used
         */
        int n_threads() const { return n_threads_; }

        /**
         * @brief Enable/disable calculation caching
         * @param enable Whether to cache fragment calculations
         */
        void set_caching(bool enable) { caching_enabled_ = enable; }

        /**
         * @brief Check if caching is enabled
         * @return True if caching enabled
         */
        bool caching_enabled() const { return caching_enabled_; }

        // State and diagnostics
        /**
         * @brief Check if calculator is ready for calculations
         * @return True if ready
         */
        bool is_ready() const;

        /**
         * @brief Get calculator status information
         * @return Status information map
         */
        std::unordered_map<std::string, std::string> get_status() const;

        /**
         * @brief Get performance statistics
         * @return Performance statistics map
         */
        std::unordered_map<std::string, double> get_performance_stats() const;

    private:
        // Configuration and interfaces
        GlobalSCFConfig config_;
        QMInterfacePtr qm_interface_;
        GlobalSCFFragmentGenerator fragment_generator_;

        // Parallelization and caching
        int n_threads_ = 0;  // 0 = auto-detect
        bool caching_enabled_ = true;

        // Callbacks
        ProgressCallback progress_callback_;
        LogCallback log_callback_;
        ConvergenceCallback convergence_callback_;

        // State management
        std::unordered_map<std::string, std::future<GlobalSCFResults>> ongoing_calculations_;
        std::unordered_map<std::string, FragmentWaveFunction> wavefunction_cache_;

        // Performance tracking
        mutable std::unordered_map<std::string, double> performance_stats_;

        // DIIS acceleration
        std::unordered_map<std::size_t, std::unique_ptr<DIISAccelerator>> diis_accelerators_;

        // Internal calculation methods
        GlobalSCFResults perform_calculation(const Molecule& molecule);

        // Fragment calculation methods
        std::vector<GlobalSCFFragment> setup_fragments(const Molecule& molecule);
        std::vector<FragmentWaveFunction> calculate_initial_fragments(
            const std::vector<GlobalSCFFragment>& fragments, const Molecule& molecule);
        
        FragmentWaveFunction calculate_fragment_iteration(
            const GlobalSCFFragment& fragment,
            const Molecule& molecule,
            const std::vector<FragmentWaveFunction>& current_wavefunctions,
            int iteration);

        // Embedding potential methods
        xt::xarray<double> generate_embedding_potential(
            const GlobalSCFFragment& target_fragment,
            const std::vector<FragmentWaveFunction>& environment_wavefunctions);

        std::vector<double> generate_point_charges(
            const GlobalSCFFragment& target_fragment,
            const std::vector<FragmentWaveFunction>& environment_wavefunctions);

        std::vector<std::array<double, 3>> get_point_charge_positions(
            const GlobalSCFFragment& target_fragment,
            const std::vector<FragmentWaveFunction>& environment_wavefunctions);

        // Polarization methods
        void update_polarization_effects(
            std::vector<FragmentWaveFunction>& wavefunctions,
            const std::vector<GlobalSCFFragment>& fragments);

        xt::xarray<double> calculate_induced_polarization(
            const GlobalSCFFragment& fragment,
            const std::vector<FragmentWaveFunction>& environment_wavefunctions);

        // Convergence checking
        bool check_convergence(const GlobalSCFIteration& iteration, const GlobalSCFConfig& config);
        GlobalSCFIteration create_iteration_data(
            int iteration_number,
            const std::vector<FragmentWaveFunction>& current_wavefunctions,
            const std::vector<FragmentWaveFunction>& previous_wavefunctions);

        double calculate_energy_change(
            const std::vector<FragmentWaveFunction>& current_wavefunctions,
            const std::vector<FragmentWaveFunction>& previous_wavefunctions);

        double calculate_density_change(
            const std::vector<FragmentWaveFunction>& current_wavefunctions,
            const std::vector<FragmentWaveFunction>& previous_wavefunctions);

        double calculate_orbital_change(
            const std::vector<FragmentWaveFunction>& current_wavefunctions,
            const std::vector<FragmentWaveFunction>& previous_wavefunctions);

        // DIIS acceleration methods
        void setup_diis_accelerators(const std::vector<GlobalSCFFragment>& fragments);
        void update_diis(int iteration, const std::vector<FragmentWaveFunction>& wavefunctions);
        void apply_diis_acceleration(std::vector<FragmentWaveFunction>& wavefunctions);

        // Utility methods
        void update_progress(int iteration, double progress, const std::string& status);
        void log_message(const std::string& message);
        void notify_convergence(const GlobalSCFIteration& iteration);

        std::string generate_cache_key(const Fragment& fragment, 
            const std::string& method, const std::string& basis_set,
            const xt::xarray<double>& embedding_potential);

        bool is_calculation_cached(const std::string& cache_key);
        FragmentWaveFunction get_cached_wavefunction(const std::string& cache_key);
        void cache_wavefunction(const std::string& cache_key, const FragmentWaveFunction& wf);

        // Validation helpers
        std::vector<std::string> validate_molecule(const Molecule& molecule);
        std::vector<std::string> validate_config(const GlobalSCFConfig& config);
        std::vector<std::string> validate_qm_interface();
        std::vector<std::string> validate_fragments(
            const std::vector<GlobalSCFFragment>& fragments, const Molecule& molecule);

        // Error handling
        void handle_calculation_error(const std::exception& e, 
            const GlobalSCFFragment& fragment, int iteration);

        // Performance tracking
        void update_performance_stats(const std::string& key, double value);
        void record_timing(const std::string& operation, 
                          std::chrono::duration<double> duration);
    };

    /**
     * @brief Factory function for creating Global SCF calculators
     * @param qm_software QM software name ("pyscf", "gaussian", "orca", etc.)
     * @param options Software-specific options
     * @return Global SCF calculator instance
     */
    std::unique_ptr<GlobalSCFCalculator> create_global_scf_calculator(
        const std::string& qm_software,
        const std::unordered_map<std::string, std::string>& options = {});

    /**
     * @brief Convenience function for simple Global SCF calculations
     * @param molecule Input molecule
     * @param embedding_type Embedding potential type
     * @param qm_method QM method string
     * @param basis_set Basis set string
     * @return Global SCF calculation results
     */
    GlobalSCFResults calculate_global_scf(const Molecule& molecule, 
        GlobalSCFConfig::EmbeddingType embedding_type = GlobalSCFConfig::EmbeddingType::COULOMB,
        const std::string& qm_method = "HF", const std::string& basis_set = "6-31G*");

}

#endif // LIBFRAG_GLOBAL_SCF_CALCULATOR_HPP
