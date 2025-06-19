#pragma once
#ifndef LIBFRAG_MBE_CONFIG_HPP
#define LIBFRAG_MBE_CONFIG_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <optional>

namespace libfrag {

    /**
     * @brief Configuration settings for Many-Body Expansion calculations
     * 
     * This class stores all parameters needed to configure MBE calculations,
     * including the expansion order, fragmentation scheme, convergence criteria,
     * and quantum chemistry method specifications.
     */
    class MBEConfig {
    public:
        /**
         * @brief Fragmentation schemes available for MBE
         */
        enum class FragmentationScheme {
            ATOMIC,           ///< Each atom is a fragment
            MOLECULAR,        ///< Each molecule is a fragment  
            FUNCTIONAL_GROUP, ///< Functional groups are fragments
            CUSTOM,           ///< User-defined fragmentation
            DISTANCE_BASED    ///< Distance-based fragmentation
        };

        /**
         * @brief Interaction truncation methods
         */
        enum class TruncationMethod {
            ORDER_BASED,      ///< Truncate at specific N-body order
            DISTANCE_BASED,   ///< Truncate based on fragment distances
            ENERGY_BASED      ///< Truncate based on energy contributions
        };

        // Constructors
        MBEConfig() = default;
        
        /**
         * @brief Construct MBE configuration with basic settings
         * @param max_order Maximum N-body order (e.g., 2 for 2-body, 3 for 3-body)
         * @param scheme Fragmentation scheme to use
         */
        MBEConfig(int max_order, FragmentationScheme scheme = FragmentationScheme::MOLECULAR);

        // Configuration setters
        /**
         * @brief Set the maximum N-body expansion order
         * @param order Maximum order (typically 2-4)
         */
        void set_max_order(int order);

        /**
         * @brief Set the fragmentation scheme
         * @param scheme Fragmentation method to use
         */
        void set_fragmentation_scheme(FragmentationScheme scheme);

        /**
         * @brief Set distance cutoff for interactions
         * @param cutoff Maximum distance between fragments (Angstroms)
         */
        void set_distance_cutoff(double cutoff);

        /**
         * @brief Set energy threshold for truncation
         * @param threshold Energy threshold in hartrees
         */
        void set_energy_threshold(double threshold);

        /**
         * @brief Set quantum chemistry method
         * @param method QM method string (e.g., "HF/6-31G*", "B3LYP/def2-TZVP")
         */
        void set_qm_method(const std::string& method);

        /**
         * @brief Set basis set
         * @param basis Basis set name
         */
        void set_basis_set(const std::string& basis);

        /**
         * @brief Enable/disable charge embedding
         * @param enable Whether to use charge embedding
         */
        void set_charge_embedding(bool enable);

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

        // Configuration getters
        int max_order() const { return max_order_; }
        FragmentationScheme fragmentation_scheme() const { return fragmentation_scheme_; }
        double distance_cutoff() const { return distance_cutoff_; }
        double energy_threshold() const { return energy_threshold_; }
        const std::string& qm_method() const { return qm_method_; }
        const std::string& basis_set() const { return basis_set_; }
        bool charge_embedding() const { return charge_embedding_; }
        TruncationMethod truncation_method() const { return truncation_method_; }
        
        const std::vector<std::vector<std::size_t>>& custom_fragments() const { 
            return custom_fragments_; 
        }
        
        const std::unordered_map<std::string, std::string>& qm_options() const { 
            return qm_options_; 
        }

        // Validation and utilities
        /**
         * @brief Validate configuration settings
         * @return True if configuration is valid
         */
        bool is_valid() const;

        /**
         * @brief Get configuration as string for logging
         * @return String representation of configuration
         */
        std::string to_string() const;

        /**
         * @brief Create configuration from JSON-like string
         * @param config_string JSON configuration string
         * @return MBEConfig object
         */
        static MBEConfig from_string(const std::string& config_string);

        // Default configurations
        /**
         * @brief Create default 2-body MBE configuration
         * @return Default 2-body config
         */
        static MBEConfig default_2body();

        /**
         * @brief Create default 3-body MBE configuration
         * @return Default 3-body config
         */
        static MBEConfig default_3body();

    private:
        // Core MBE settings
        int max_order_ = 2;
        FragmentationScheme fragmentation_scheme_ = FragmentationScheme::MOLECULAR;
        TruncationMethod truncation_method_ = TruncationMethod::ORDER_BASED;
        
        // Cutoff criteria
        double distance_cutoff_ = 10.0;  // Angstroms
        double energy_threshold_ = 1e-6; // Hartrees
        
        // Quantum chemistry settings
        std::string qm_method_ = "HF";
        std::string basis_set_ = "6-31G*";
        bool charge_embedding_ = false;
        
        // Custom fragmentation
        std::vector<std::vector<std::size_t>> custom_fragments_;
        
        // Additional QM options
        std::unordered_map<std::string, std::string> qm_options_;
        
        // Validation helpers
        bool validate_order() const;
        bool validate_cutoffs() const;
        bool validate_qm_settings() const;
    };

    /**
     * @brief Convert fragmentation scheme to string
     */
    std::string to_string(MBEConfig::FragmentationScheme scheme);

    /**
     * @brief Convert truncation method to string  
     */
    std::string to_string(MBEConfig::TruncationMethod method);

} // namespace libfrag

#endif // LIBFRAG_MBE_CONFIG_HPP
