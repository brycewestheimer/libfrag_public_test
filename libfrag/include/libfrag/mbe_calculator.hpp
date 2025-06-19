#pragma once
#ifndef LIBFRAG_MBE_CALCULATOR_HPP
#define LIBFRAG_MBE_CALCULATOR_HPP

#include "libfrag/molecule.hpp"
#include "libfrag/fragment.hpp"
#include "libfrag/mbe_config.hpp"
#include "libfrag/mbe_results.hpp"
#include "libfrag/mbe_fragment_generator.hpp"
#include <memory>
#include <functional>
#include <future>
#include <string>

namespace libfrag {

    /**
     * @brief Quantum chemistry calculation interface for MBE
     * 
     * Abstract base class defining interface for quantum chemistry calculations
     * needed by the MBE calculator. Users can implement this interface to connect
     * with their preferred quantum chemistry software (PySCF, Gaussian, ORCA, etc.)
     */
    class QMCalculatorInterface {
    public:
        virtual ~QMCalculatorInterface() = default;

        /**
         * @brief Perform single-point energy calculation
         * @param fragment Fragment to calculate
         * @param method QM method string
         * @param basis_set Basis set string
         * @param charge Total charge
         * @param multiplicity Spin multiplicity
         * @param options Additional QM options
         * @return Fragment calculation result
         */
        virtual FragmentCalculationResult calculate_energy(
            const Fragment& fragment,
            const std::string& method,
            const std::string& basis_set,
            int charge = 0,
            int multiplicity = 1,
            const std::unordered_map<std::string, std::string>& options = {}) = 0;

        /**
         * @brief Check if QM software is available
         * @return True if QM software is accessible
         */
        virtual bool is_available() const = 0;

        /**
         * @brief Get QM software name/version
         * @return Software identification string
         */
        virtual std::string software_info() const = 0;
    };

    /**
     * @brief Main Many-Body Expansion calculator
     * 
     * This class orchestrates the complete MBE calculation workflow:
     * 1. Fragment generation based on configuration
     * 2. Parallel quantum chemistry calculations for all fragment combinations
     * 3. Energy aggregation and error analysis
     * 4. Result compilation and export
     */
    class MBECalculator {
    public:
        // Type aliases
        using QMCalculatorPtr = std::unique_ptr<QMCalculatorInterface>;
        using ProgressCallback = std::function<void(double progress, const std::string& status)>;
        using LogCallback = std::function<void(const std::string& message)>;

        // Constructors
        MBECalculator() = default;
        
        /**
         * @brief Construct MBE calculator with QM interface
         * @param qm_calculator Quantum chemistry calculator implementation
         */
        explicit MBECalculator(QMCalculatorPtr qm_calculator);

        /**
         * @brief Construct MBE calculator with QM interface and configuration
         * @param qm_calculator Quantum chemistry calculator implementation
         * @param config MBE configuration
         */
        MBECalculator(QMCalculatorPtr qm_calculator, const MBEConfig& config);

        // Configuration
        /**
         * @brief Set MBE configuration
         * @param config MBE configuration
         */
        void set_config(const MBEConfig& config) { config_ = config; }

        /**
         * @brief Get current MBE configuration
         * @return Current configuration
         */
        const MBEConfig& config() const { return config_; }

        /**
         * @brief Set quantum chemistry calculator
         * @param qm_calculator QM calculator implementation
         */
        void set_qm_calculator(QMCalculatorPtr qm_calculator);

        // Callback setup
        /**
         * @brief Set progress callback for monitoring calculations
         * @param callback Progress callback function
         */
        void set_progress_callback(ProgressCallback callback) { progress_callback_ = callback; }

        /**
         * @brief Set logging callback for calculation messages
         * @param callback Log callback function
         */
        void set_log_callback(LogCallback callback) { log_callback_ = callback; }

        // Main calculation methods
        /**
         * @brief Perform complete MBE calculation
         * @param molecule Input molecule
         * @return MBE calculation results
         */
        MBEResults calculate(const Molecule& molecule);

        /**
         * @brief Perform MBE calculation with custom configuration
         * @param molecule Input molecule
         * @param config MBE configuration
         * @return MBE calculation results
         */
        MBEResults calculate(const Molecule& molecule, const MBEConfig& config);

        /**
         * @brief Perform async MBE calculation
         * @param molecule Input molecule
         * @return Future for MBE calculation results
         */
        std::future<MBEResults> calculate_async(const Molecule& molecule);

        // Incremental calculation methods
        /**
         * @brief Start incremental MBE calculation
         * @param molecule Input molecule
         * @return Calculation handle for incremental updates
         */
        std::string start_calculation(const Molecule& molecule);

        /**
         * @brief Get partial results from ongoing calculation
         * @param calculation_id Calculation handle
         * @return Partial MBE results
         */
        MBEResults get_partial_results(const std::string& calculation_id);

        /**
         * @brief Check if calculation is complete
         * @param calculation_id Calculation handle
         * @return True if calculation finished
         */
        bool is_calculation_complete(const std::string& calculation_id);

        /**
         * @brief Cancel ongoing calculation
         * @param calculation_id Calculation handle
         */
        void cancel_calculation(const std::string& calculation_id);

        // Analysis methods
        /**
         * @brief Estimate computational cost
         * @param molecule Input molecule
         * @param config MBE configuration
         * @return Estimated cost metrics
         */
        std::unordered_map<std::string, double> estimate_cost(
            const Molecule& molecule, const MBEConfig& config);

        /**
         * @brief Preview fragment combinations without calculation
         * @param molecule Input molecule
         * @param config MBE configuration
         * @return Fragment combination preview
         */
        std::unordered_map<int, std::size_t> preview_fragments(
            const Molecule& molecule, const MBEConfig& config);

        /**
         * @brief Validate calculation setup
         * @param molecule Input molecule
         * @param config MBE configuration
         * @return Validation report
         */
        std::vector<std::string> validate_setup(
            const Molecule& molecule, const MBEConfig& config);

        // Parallelization control
        /**
         * @brief Set number of parallel threads
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
         * @brief Clear calculation cache
         */
        void clear_cache();

        /**
         * @brief Get performance statistics
         * @return Performance metrics
         */
        std::unordered_map<std::string, double> get_performance_stats() const;

    private:
        // Core components
        MBEConfig config_;
        QMCalculatorPtr qm_calculator_;
        MBEFragmentGenerator fragment_generator_;
        
        // Parallelization
        int n_threads_ = 0;  // 0 = auto-detect
        bool caching_enabled_ = true;
        
        // Callbacks
        ProgressCallback progress_callback_;
        LogCallback log_callback_;
        
        // State management
        std::unordered_map<std::string, std::future<MBEResults>> ongoing_calculations_;
        std::unordered_map<std::string, FragmentCalculationResult> calculation_cache_;
        
        // Performance tracking
        mutable std::unordered_map<std::string, double> performance_stats_;
        
        // Internal calculation methods
        MBEResults perform_calculation(const Molecule& molecule);
        
        std::vector<FragmentCalculationResult> calculate_n_body_terms(
            const std::vector<std::size_t>& fragment_indices,
            const std::vector<std::shared_ptr<Fragment>>& fragments,
            int n_body);
        
        FragmentCalculationResult calculate_single_fragment(
            const std::vector<std::size_t>& fragment_indices,
            const std::vector<std::shared_ptr<Fragment>>& fragments);
        
        // Utility methods
        void update_progress(double progress, const std::string& status);
        void log_message(const std::string& message);
        
        std::string generate_cache_key(const Fragment& fragment, 
            const std::string& method, const std::string& basis_set,
            int charge, int multiplicity);
        
        bool is_calculation_cached(const std::string& cache_key);
        FragmentCalculationResult get_cached_result(const std::string& cache_key);
        void cache_result(const std::string& cache_key, const FragmentCalculationResult& result);
        
        // Validation helpers
        std::vector<std::string> validate_molecule(const Molecule& molecule);
        std::vector<std::string> validate_config(const MBEConfig& config);
        std::vector<std::string> validate_qm_calculator();
        
        // Error handling
        void handle_calculation_error(const std::exception& e, 
            const std::vector<std::size_t>& fragment_indices);
    };

    /**
     * @brief Factory function for creating MBE calculators
     * @param qm_software QM software name ("pyscf", "gaussian", "orca", etc.)
     * @param options Software-specific options
     * @return MBE calculator instance
     */
    std::unique_ptr<MBECalculator> create_mbe_calculator(
        const std::string& qm_software,
        const std::unordered_map<std::string, std::string>& options = {});

    /**
     * @brief Convenience function for simple MBE calculations
     * @param molecule Input molecule
     * @param max_order Maximum N-body order
     * @param qm_method QM method string
     * @param basis_set Basis set string
     * @return MBE calculation results
     */
    MBEResults calculate_mbe(const Molecule& molecule, int max_order,
        const std::string& qm_method = "HF", const std::string& basis_set = "6-31G*");

} // namespace libfrag

#endif // LIBFRAG_MBE_CALCULATOR_HPP
