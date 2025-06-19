#pragma once
#ifndef LIBFRAG_GLOBAL_SCF_CONFIG_HPP
#define LIBFRAG_GLOBAL_SCF_CONFIG_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

namespace libfrag {

    /**
     * @brief Configuration settings for Global SCF calculations
     * 
     * This class stores all parameters needed to configure Global SCF calculations,
     * including the fragmentation scheme, embedding potential types, convergence criteria,
     * and quantum chemistry method specifications for fragment-based SCF methods like FMO.
     */
    class GlobalSCFConfig {
    public:
        /**
         * @brief Fragmentation schemes available for Global SCF
         */
        enum class FragmentationScheme {
            ATOMIC,           ///< Each atom is a fragment
            MOLECULAR,        ///< Each molecule is a fragment  
            FUNCTIONAL_GROUP, ///< Functional groups are fragments
            CUSTOM,           ///< User-defined fragmentation
            DISTANCE_BASED,   ///< Distance-based fragmentation
            FMO_LIKE         ///< Fragment Molecular Orbital style fragmentation
        };

        /**
         * @brief Types of embedding potentials
         */
        enum class EmbeddingType {
            NONE,             ///< No embedding (in vacuo calculations)
            COULOMB,          ///< Coulomb embedding (point charges)
            POLARIZABLE,      ///< Polarizable embedding
            DENSITY_BASED,    ///< Density-based embedding
            COMBINED          ///< Combined Coulomb + polarizable
        };

        /**
         * @brief SCF convergence criteria
         */
        enum class ConvergenceType {
            ENERGY,           ///< Energy-based convergence
            DENSITY,          ///< Density matrix convergence
            ORBITAL,          ///< Orbital coefficient convergence
            COMBINED          ///< Combined criteria
        };

        // Constructors
        GlobalSCFConfig() = default;
        
        /**
         * @brief Construct with basic settings
         * @param scheme Fragmentation scheme
         * @param embedding Embedding potential type
         */
        GlobalSCFConfig(FragmentationScheme scheme, EmbeddingType embedding = EmbeddingType::COULOMB);

        // Configuration setters
        /**
         * @brief Set fragmentation scheme
         * @param scheme Fragmentation method
         */
        void set_fragmentation_scheme(FragmentationScheme scheme);

        /**
         * @brief Set embedding potential type
         * @param embedding Type of embedding potential
         */
        void set_embedding_type(EmbeddingType embedding);

        /**
         * @brief Set convergence criteria type
         * @param convergence Type of convergence criteria
         */
        void set_convergence_type(ConvergenceType convergence);

        /**
         * @brief Set maximum SCF iterations
         * @param max_iter Maximum number of global SCF iterations
         */
        void set_max_scf_iterations(int max_iter);

        /**
         * @brief Set energy convergence threshold
         * @param threshold Energy convergence threshold in hartrees
         */
        void set_energy_threshold(double threshold);

        /**
         * @brief Set density convergence threshold
         * @param threshold Density matrix convergence threshold
         */
        void set_density_threshold(double threshold);

        /**
         * @brief Set orbital convergence threshold
         * @param threshold Orbital coefficient convergence threshold
         */
        void set_orbital_threshold(double threshold);

        /**
         * @brief Set quantum chemistry method
         * @param method QM method string (e.g., "HF", "B3LYP", "MP2")
         */
        void set_qm_method(const std::string& method);

        /**
         * @brief Set basis set
         * @param basis Basis set name
         */
        void set_basis_set(const std::string& basis);

        /**
         * @brief Enable/disable DIIS acceleration
         * @param enable Whether to use DIIS for SCF acceleration
         */
        void set_diis_enabled(bool enable);

        /**
         * @brief Set DIIS subspace size
         * @param size Maximum number of DIIS vectors
         */
        void set_diis_size(int size);

        /**
         * @brief Set fragment buffer distance
         * @param distance Buffer distance around fragments in Angstroms
         */
        void set_buffer_distance(double distance);

        /**
         * @brief Set custom fragment definitions
         * @param fragments Vector of atom indices for each fragment
         */
        void set_custom_fragments(const std::vector<std::vector<std::size_t>>& fragments);

        /**
         * @brief Set additional QM software options
         * @param options Key-value pairs for QM software
         */
        void set_qm_options(const std::unordered_map<std::string, std::string>& options);

        /**
         * @brief Enable/disable fragment polarization
         * @param enable Whether fragments can polarize each other
         */
        void set_polarization_enabled(bool enable);

        /**
         * @brief Set polarization convergence threshold
         * @param threshold Convergence threshold for polarization
         */
        void set_polarization_threshold(double threshold);

        // Configuration getters
        FragmentationScheme fragmentation_scheme() const { return fragmentation_scheme_; }
        EmbeddingType embedding_type() const { return embedding_type_; }
        ConvergenceType convergence_type() const { return convergence_type_; }
        int max_scf_iterations() const { return max_scf_iterations_; }
        double energy_threshold() const { return energy_threshold_; }
        double density_threshold() const { return density_threshold_; }
        double orbital_threshold() const { return orbital_threshold_; }
        const std::string& qm_method() const { return qm_method_; }
        const std::string& basis_set() const { return basis_set_; }
        bool diis_enabled() const { return diis_enabled_; }
        int diis_size() const { return diis_size_; }
        double buffer_distance() const { return buffer_distance_; }
        bool polarization_enabled() const { return polarization_enabled_; }
        double polarization_threshold() const { return polarization_threshold_; }
        
        const std::vector<std::vector<std::size_t>>& custom_fragments() const { 
            return custom_fragments_; 
        }
        const std::unordered_map<std::string, std::string>& qm_options() const { 
            return qm_options_; 
        }

        // Validation and factory methods
        /**
         * @brief Validate configuration settings
         * @return True if configuration is valid
         */
        bool is_valid() const;

        /**
         * @brief Get detailed configuration string
         * @return Formatted configuration description
         */
        std::string to_string() const;

        /**
         * @brief Create FMO-style configuration
         * @param method QM method
         * @param basis Basis set
         * @return FMO configuration
         */
        static GlobalSCFConfig create_fmo_config(const std::string& method = "HF", 
                                                  const std::string& basis = "6-31G*");

        /**
         * @brief Create simple Coulomb embedding configuration
         * @param method QM method
         * @param basis Basis set
         * @return Coulomb embedding configuration
         */
        static GlobalSCFConfig create_coulomb_config(const std::string& method = "HF", 
                                                      const std::string& basis = "6-31G*");

        /**
         * @brief Create polarizable embedding configuration
         * @param method QM method
         * @param basis Basis set
         * @return Polarizable embedding configuration
         */
        static GlobalSCFConfig create_polarizable_config(const std::string& method = "HF", 
                                                          const std::string& basis = "6-31G*");

        // String conversion utilities
        static std::string fragmentation_scheme_to_string(FragmentationScheme scheme);
        static std::string embedding_type_to_string(EmbeddingType embedding);
        static std::string convergence_type_to_string(ConvergenceType convergence);

    private:
        // Core configuration
        FragmentationScheme fragmentation_scheme_ = FragmentationScheme::MOLECULAR;
        EmbeddingType embedding_type_ = EmbeddingType::COULOMB;
        ConvergenceType convergence_type_ = ConvergenceType::COMBINED;
        
        // SCF parameters
        int max_scf_iterations_ = 50;
        double energy_threshold_ = 1e-6; // Hartrees
        double density_threshold_ = 1e-5; // Density matrix convergence
        double orbital_threshold_ = 1e-5; // Orbital coefficient convergence
        
        // DIIS parameters
        bool diis_enabled_ = true;
        int diis_size_ = 8;
        
        // Polarization parameters
        bool polarization_enabled_ = true;
        double polarization_threshold_ = 1e-6;
        
        // Fragment parameters
        double buffer_distance_ = 2.0; // Angstroms
        
        // Quantum chemistry settings
        std::string qm_method_ = "HF";
        std::string basis_set_ = "6-31G*";
        
        // Custom fragmentation
        std::vector<std::vector<std::size_t>> custom_fragments_;
        
        // Additional QM options
        std::unordered_map<std::string, std::string> qm_options_;
        
        // Validation helpers
        bool validate_thresholds() const;
        bool validate_qm_settings() const;
        bool validate_diis_settings() const;
    };

}

#endif // LIBFRAG_GLOBAL_SCF_CONFIG_HPP
