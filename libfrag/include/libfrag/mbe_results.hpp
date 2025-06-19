#pragma once
#ifndef LIBFRAG_MBE_RESULTS_HPP
#define LIBFRAG_MBE_RESULTS_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <chrono>
#include "xtensor/xarray.hpp"

namespace libfrag {

    // Forward declarations
    class Fragment;
    class MBEConfig;

    /**
     * @brief Individual fragment calculation result
     * 
     * Stores the quantum chemistry calculation results for a single fragment
     * or fragment combination in the MBE expansion.
     */
    struct FragmentCalculationResult {
        // Fragment identification
        std::vector<std::size_t> fragment_indices;  ///< Indices of fragments involved
        int n_body_order;                           ///< N-body order (1, 2, 3, etc.)
        std::string fragment_id;                    ///< Unique identifier
        
        // Energies (all in hartrees)
        double total_energy = 0.0;                  ///< Total electronic energy
        double nuclear_repulsion = 0.0;             ///< Nuclear repulsion energy
        double electronic_energy = 0.0;             ///< Electronic energy
        double correlation_energy = 0.0;            ///< Correlation energy (if applicable)
        
        // Properties
        std::unordered_map<std::string, double> properties;  ///< Additional properties
        
        // Quantum chemistry details
        std::string qm_method;                      ///< QM method used
        std::string basis_set;                      ///< Basis set used
        int n_electrons = 0;                        ///< Number of electrons
        int multiplicity = 1;                       ///< Spin multiplicity
        int charge = 0;                             ///< Total charge
        
        // Computational details
        std::chrono::duration<double> computation_time; ///< Wall clock time
        bool converged = false;                     ///< SCF convergence status
        int scf_iterations = 0;                     ///< Number of SCF iterations
        std::string status;                         ///< Calculation status message
        
        // Optional wavefunction data
        std::optional<xt::xarray<double>> density_matrix;    ///< Density matrix
        std::optional<xt::xarray<double>> mo_coefficients;   ///< Molecular orbital coefficients
        std::optional<xt::xarray<double>> mo_energies;       ///< Molecular orbital energies
        
        // Constructor
        FragmentCalculationResult() = default;
        FragmentCalculationResult(const std::vector<std::size_t>& indices, int order);
        
        // Utilities
        std::string to_string() const;
        bool is_valid() const;
    };

    /**
     * @brief Complete Many-Body Expansion calculation results
     * 
     * Stores all results from an MBE calculation including individual fragment
     * energies, interaction energies, and the final total energy with error analysis.
     */
    class MBEResults {
    public:
        // Constructors
        MBEResults() = default;
        MBEResults(const MBEConfig& config);

        // Result storage
        /**
         * @brief Add a fragment calculation result
         * @param result Individual fragment calculation result
         */
        void add_fragment_result(const FragmentCalculationResult& result);

        /**
         * @brief Add fragment result with move semantics
         * @param result Fragment calculation result (moved)
         */
        void add_fragment_result(FragmentCalculationResult&& result);

        // Energy accessors
        /**
         * @brief Get total MBE energy
         * @return Total energy in hartrees
         */
        double total_energy() const { return total_energy_; }

        /**
         * @brief Get energy contributions by N-body order
         * @param order N-body order (1, 2, 3, etc.)
         * @return Energy contribution for that order
         */
        double energy_contribution(int order) const;

        /**
         * @brief Get all energy contributions by order
         * @return Map of order -> energy contribution
         */
        std::unordered_map<int, double> energy_contributions() const;

        /**
         * @brief Get individual fragment energies (1-body terms)
         * @return Vector of 1-body energies
         */
        std::vector<double> fragment_energies() const;

        /**
         * @brief Get interaction energies by order
         * @param order N-body interaction order
         * @return Vector of interaction energies
         */
        std::vector<double> interaction_energies(int order) const;

        // Result accessors
        /**
         * @brief Get all fragment calculation results
         * @return Vector of all fragment results
         */
        const std::vector<FragmentCalculationResult>& fragment_results() const { 
            return fragment_results_; 
        }

        /**
         * @brief Get fragment results by N-body order
         * @param order N-body order
         * @return Vector of results for that order
         */
        std::vector<FragmentCalculationResult> results_by_order(int order) const;

        /**
         * @brief Get maximum N-body order calculated
         * @return Maximum order
         */
        int max_order() const { return max_order_; }

        /**
         * @brief Get number of fragments
         * @return Number of fragments
         */
        std::size_t n_fragments() const { return n_fragments_; }

        // Error analysis
        /**
         * @brief Estimate MBE truncation error
         * @return Estimated error in hartrees
         */
        double estimated_truncation_error() const;

        /**
         * @brief Check if MBE series has converged
         * @param threshold Convergence threshold
         * @return True if converged
         */
        bool is_converged(double threshold = 1e-6) const;

        /**
         * @brief Get convergence analysis
         * @return Map with convergence metrics
         */
        std::unordered_map<std::string, double> convergence_analysis() const;

        // Computational statistics
        /**
         * @brief Get total computation time
         * @return Total wall clock time
         */
        std::chrono::duration<double> total_computation_time() const;

        /**
         * @brief Get computation time by N-body order
         * @param order N-body order
         * @return Computation time for that order
         */
        std::chrono::duration<double> computation_time_by_order(int order) const;

        /**
         * @brief Get performance statistics
         * @return Map with performance metrics
         */
        std::unordered_map<std::string, double> performance_statistics() const;

        // Export and analysis
        /**
         * @brief Export results to JSON string
         * @return JSON representation of results
         */
        std::string to_json() const;

        /**
         * @brief Export results to CSV format
         * @return CSV representation of results
         */
        std::string to_csv() const;

        /**
         * @brief Create detailed summary report
         * @return Formatted summary string
         */
        std::string summary_report() const;

        /**
         * @brief Export energy decomposition table
         * @return Formatted table of energy contributions
         */
        std::string energy_decomposition_table() const;

        // Validation and utilities
        /**
         * @brief Validate results consistency
         * @return True if results are consistent
         */
        bool validate() const;

        /**
         * @brief Clear all results
         */
        void clear();

        /**
         * @brief Check if results are empty
         * @return True if no results stored
         */
        bool empty() const { return fragment_results_.empty(); }

        /**
         * @brief Get number of calculations performed
         * @return Number of fragment calculations
         */
        std::size_t n_calculations() const { return fragment_results_.size(); }

        // Operators
        /**
         * @brief Access fragment result by index
         * @param index Result index
         * @return Reference to fragment result
         */
        const FragmentCalculationResult& operator[](std::size_t index) const;

    private:
        // Core results data
        std::vector<FragmentCalculationResult> fragment_results_;
        
        // Computed totals
        double total_energy_ = 0.0;
        int max_order_ = 0;
        std::size_t n_fragments_ = 0;
        
        // Energy contributions by order
        std::unordered_map<int, double> energy_by_order_;
        
        // Configuration used
        std::shared_ptr<MBEConfig> config_;
        
        // Computation metadata
        std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
        std::chrono::time_point<std::chrono::high_resolution_clock> end_time_;
        
        // Internal methods
        void update_totals();
        void compute_energy_contributions();
        double calculate_interaction_energy(const std::vector<std::size_t>& fragment_indices) const;
        std::string format_energy_table() const;
        std::string format_timing_table() const;
    };

    /**
     * @brief Comparison operator for fragment calculation results
     */
    bool operator==(const FragmentCalculationResult& lhs, const FragmentCalculationResult& rhs);

    /**
     * @brief Stream output operator for fragment results
     */
    std::ostream& operator<<(std::ostream& os, const FragmentCalculationResult& result);

    /**
     * @brief Stream output operator for MBE results
     */
    std::ostream& operator<<(std::ostream& os, const MBEResults& results);

} // namespace libfrag

#endif // LIBFRAG_MBE_RESULTS_HPP
