#pragma once
#ifndef LIBFRAG_GLOBAL_SCF_RESULTS_HPP
#define LIBFRAG_GLOBAL_SCF_RESULTS_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <chrono>
#include <memory>
#include <optional>
#include "xtensor/xarray.hpp"
#include "global_scf_config.hpp"

namespace libfrag {

    /**
     * @brief Fragment wave function data
     * 
     * Stores the quantum chemistry calculation results and wave function data
     * for a single fragment in the Global SCF calculation.
     */
    struct FragmentWaveFunction {
        // Fragment identification
        std::size_t fragment_id;                    ///< Fragment identifier
        std::vector<std::size_t> atom_indices;     ///< Atom indices in fragment
        std::string fragment_name;                  ///< Human-readable fragment name
        
        // Energies (all in hartrees)
        double total_energy = 0.0;                 ///< Total fragment energy
        double electronic_energy = 0.0;            ///< Electronic energy
        double nuclear_repulsion = 0.0;            ///< Nuclear repulsion energy
        double embedding_energy = 0.0;             ///< Energy from embedding potential
        double polarization_energy = 0.0;          ///< Polarization energy
        
        // Wave function data
        xt::xarray<double> density_matrix;         ///< Density matrix
        xt::xarray<double> mo_coefficients;        ///< Molecular orbital coefficients
        xt::xarray<double> mo_energies;            ///< Molecular orbital energies
        xt::xarray<double> mo_occupancies;         ///< Molecular orbital occupancies
        
        // Embedding potential data
        xt::xarray<double> embedding_potential;    ///< Embedding potential matrix
        std::vector<double> point_charges;         ///< Point charges from environment
        std::vector<std::array<double, 3>> charge_positions; ///< Positions of point charges
        
        // Polarization data
        xt::xarray<double> polarization_potential; ///< Polarization potential
        xt::xarray<double> induced_dipoles;        ///< Induced dipole moments
        
        // Computational details
        std::string qm_method;                     ///< QM method used
        std::string basis_set;                     ///< Basis set used
        int n_electrons = 0;                       ///< Number of electrons
        int multiplicity = 1;                      ///< Spin multiplicity
        int charge = 0;                            ///< Total charge
        
        // Convergence information
        bool converged = false;                    ///< SCF convergence status
        int scf_iterations = 0;                    ///< Number of SCF iterations
        double final_energy_change = 0.0;         ///< Final energy change
        double final_density_change = 0.0;        ///< Final density change
        std::chrono::duration<double> computation_time; ///< Wall clock time
        
        // Status and diagnostics
        std::string status;                        ///< Calculation status message
        std::unordered_map<std::string, double> properties; ///< Additional properties
        
        // Constructor
        FragmentWaveFunction() = default;
        FragmentWaveFunction(std::size_t id, const std::vector<std::size_t>& atoms);
        
        // Utilities
        std::string to_string() const;
        bool is_valid() const;
        double get_property(const std::string& name, double default_value = 0.0) const;
        void set_property(const std::string& name, double value);
        
        // Wave function analysis
        double electronic_charge() const;
        double nuclear_charge() const;
        std::array<double, 3> electric_dipole_moment() const;
        double total_electron_density() const;
    };

    /**
     * @brief Global SCF iteration data
     * 
     * Stores information about a single global SCF iteration,
     * including energies, convergence criteria, and timing.
     */
    struct GlobalSCFIteration {
        int iteration_number = 0;                  ///< Iteration number (1-based)
        
        // Energies for this iteration
        double total_energy = 0.0;                 ///< Total system energy
        double electronic_energy = 0.0;            ///< Total electronic energy
        double nuclear_repulsion = 0.0;            ///< Total nuclear repulsion
        double embedding_energy = 0.0;             ///< Total embedding energy
        double polarization_energy = 0.0;          ///< Total polarization energy
        double interaction_energy = 0.0;           ///< Inter-fragment interaction energy
        
        // Convergence criteria
        double energy_change = 0.0;                ///< Energy change from previous iteration
        double max_density_change = 0.0;           ///< Maximum density matrix change
        double rms_density_change = 0.0;           ///< RMS density matrix change
        double max_orbital_change = 0.0;           ///< Maximum orbital coefficient change
        double polarization_change = 0.0;          ///< Polarization energy change
        
        // Timing information
        std::chrono::duration<double> iteration_time; ///< Time for this iteration
        std::chrono::duration<double> fragment_time;  ///< Time for fragment calculations
        std::chrono::duration<double> embedding_time; ///< Time for embedding updates
        
        // Convergence status
        bool energy_converged = false;             ///< Energy convergence status
        bool density_converged = false;            ///< Density convergence status
        bool orbital_converged = false;            ///< Orbital convergence status
        bool polarization_converged = false;       ///< Polarization convergence status
        bool overall_converged = false;            ///< Overall convergence status
        
        // DIIS information
        bool diis_used = false;                    ///< Whether DIIS was used
        int diis_dimension = 0;                    ///< DIIS subspace dimension
        double diis_error = 0.0;                   ///< DIIS error norm
        
        // Diagnostics
        std::string status_message;                ///< Status message for this iteration
        std::unordered_map<std::string, double> diagnostics; ///< Additional diagnostics
        
        // Constructor
        GlobalSCFIteration() = default;
        GlobalSCFIteration(int iter_num);
        
        // Utilities
        std::string to_string() const;
        bool is_converged(const GlobalSCFConfig& config) const;
        void set_diagnostic(const std::string& name, double value);
        double get_diagnostic(const std::string& name, double default_value = 0.0) const;
    };

    /**
     * @brief Results container for Global SCF calculations
     * 
     * This class stores the complete results of a Global SCF calculation,
     * including fragment wave functions, iteration history, convergence analysis,
     * and performance statistics.
     */
    class GlobalSCFResults {
    public:
        // Constructors
        GlobalSCFResults() = default;
        explicit GlobalSCFResults(const GlobalSCFConfig& config);

        // Fragment wave function management
        /**
         * @brief Add fragment wave function
         * @param wf Fragment wave function to add
         */
        void add_fragment_wavefunction(FragmentWaveFunction wf);

        /**
         * @brief Get fragment wave function by ID
         * @param fragment_id Fragment identifier
         * @return Fragment wave function reference
         */
        const FragmentWaveFunction& get_fragment_wavefunction(std::size_t fragment_id) const;

        /**
         * @brief Get fragment wave function by ID (mutable)
         * @param fragment_id Fragment identifier
         * @return Fragment wave function reference
         */
        FragmentWaveFunction& get_fragment_wavefunction(std::size_t fragment_id);

        /**
         * @brief Get all fragment wave functions
         * @return Vector of fragment wave functions
         */
        const std::vector<FragmentWaveFunction>& fragment_wavefunctions() const { 
            return fragment_wavefunctions_; 
        }

        // Iteration history management
        /**
         * @brief Add SCF iteration data
         * @param iteration SCF iteration data
         */
        void add_scf_iteration(GlobalSCFIteration iteration);

        /**
         * @brief Get SCF iteration by number
         * @param iteration_number Iteration number (1-based)
         * @return SCF iteration data
         */
        const GlobalSCFIteration& get_scf_iteration(int iteration_number) const;

        /**
         * @brief Get all SCF iterations
         * @return Vector of SCF iterations
         */
        const std::vector<GlobalSCFIteration>& scf_iterations() const { 
            return scf_iterations_; 
        }

        // Energy analysis
        /**
         * @brief Get final total energy
         * @return Final converged energy in hartrees
         */
        double final_energy() const;

        /**
         * @brief Get energy breakdown by component
         * @return Map of energy components
         */
        std::unordered_map<std::string, double> energy_breakdown() const;

        /**
         * @brief Get fragment energies
         * @return Map of fragment ID to energy
         */
        std::unordered_map<std::size_t, double> fragment_energies() const;

        /**
         * @brief Get interaction energy between fragments
         * @return Total inter-fragment interaction energy
         */
        double interaction_energy() const;

        // Convergence analysis
        /**
         * @brief Check if calculation converged
         * @return True if calculation converged
         */
        bool is_converged() const;

        /**
         * @brief Get number of SCF iterations performed
         * @return Number of iterations
         */
        int n_iterations() const { return static_cast<int>(scf_iterations_.size()); }

        /**
         * @brief Get convergence analysis
         * @return Map of convergence metrics
         */
        std::unordered_map<std::string, double> convergence_analysis() const;

        /**
         * @brief Get convergence history for plotting
         * @param criterion Convergence criterion name
         * @return Vector of convergence values by iteration
         */
        std::vector<double> convergence_history(const std::string& criterion) const;

        // Performance analysis
        /**
         * @brief Get total computation time
         * @return Total wall clock time
         */
        std::chrono::duration<double> total_time() const;

        /**
         * @brief Get performance breakdown
         * @return Map of timing components
         */
        std::unordered_map<std::string, double> performance_breakdown() const;

        /**
         * @brief Get performance statistics
         * @return Map of performance metrics
         */
        std::unordered_map<std::string, double> performance_statistics() const;

        // Export and analysis
        /**
         * @brief Export results to JSON
         * @param filename Output filename
         * @param include_wavefunctions Whether to include wave function data
         */
        void export_json(const std::string& filename, bool include_wavefunctions = false) const;

        /**
         * @brief Export convergence data to CSV
         * @param filename Output filename
         */
        void export_convergence_csv(const std::string& filename) const;

        /**
         * @brief Export energy breakdown to CSV
         * @param filename Output filename
         */
        void export_energy_csv(const std::string& filename) const;

        /**
         * @brief Generate formatted summary report
         * @param detailed Whether to include detailed information
         * @return Formatted summary string
         */
        std::string generate_summary(bool detailed = false) const;

        // Validation and diagnostics
        /**
         * @brief Validate results consistency
         * @return Vector of validation issues (empty if valid)
         */
        std::vector<std::string> validate() const;

        /**
         * @brief Get calculation statistics
         * @return Map of calculation statistics
         */
        std::unordered_map<std::string, double> calculation_statistics() const;

        // Configuration access
        /**
         * @brief Get configuration used for calculation
         * @return Configuration object
         */
        const GlobalSCFConfig& config() const { return config_; }

        /**
         * @brief Set configuration
         * @param config Configuration object
         */
        void set_config(const GlobalSCFConfig& config) { config_ = config; }

        // Utility methods
        /**
         * @brief Get number of fragments
         * @return Number of fragments
         */
        std::size_t n_fragments() const { return fragment_wavefunctions_.size(); }

        /**
         * @brief Check if results are empty
         * @return True if no results stored
         */
        bool empty() const { return fragment_wavefunctions_.empty() && scf_iterations_.empty(); }

        /**
         * @brief Clear all results
         */
        void clear();

        /**
         * @brief Get string representation
         * @return String representation of results
         */
        std::string to_string() const;

    private:
        // Configuration
        GlobalSCFConfig config_;
        
        // Fragment data
        std::vector<FragmentWaveFunction> fragment_wavefunctions_;
        std::unordered_map<std::size_t, std::size_t> fragment_id_to_index_;
        
        // Iteration history
        std::vector<GlobalSCFIteration> scf_iterations_;
        
        // Timing information
        std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
        std::chrono::time_point<std::chrono::high_resolution_clock> end_time_;
        
        // Utility methods
        void update_fragment_index(std::size_t fragment_id, std::size_t index);
        std::size_t get_fragment_index(std::size_t fragment_id) const;
        
        // Export helpers
        std::string export_json_string(bool include_wavefunctions) const;
        void write_csv_header(std::ostream& os, const std::vector<std::string>& headers) const;
        void write_csv_row(std::ostream& os, const std::vector<double>& values) const;
    };

}

#endif // LIBFRAG_GLOBAL_SCF_RESULTS_HPP
