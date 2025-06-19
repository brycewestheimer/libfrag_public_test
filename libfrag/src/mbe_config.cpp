#include "libfrag/mbe_config.hpp"
#include <sstream>
#include <stdexcept>

namespace libfrag {

    MBEConfig::MBEConfig(int max_order, FragmentationScheme scheme) 
        : max_order_(max_order), fragmentation_scheme_(scheme) {
        if (max_order <= 0) {
            throw std::invalid_argument("Maximum order must be positive");
        }
    }

    void MBEConfig::set_max_order(int order) {
        if (order <= 0) {
            throw std::invalid_argument("Maximum order must be positive");
        }
        max_order_ = order;
    }

    void MBEConfig::set_fragmentation_scheme(FragmentationScheme scheme) {
        fragmentation_scheme_ = scheme;
    }

    void MBEConfig::set_distance_cutoff(double cutoff) {
        if (cutoff <= 0.0) {
            throw std::invalid_argument("Distance cutoff must be positive");
        }
        distance_cutoff_ = cutoff;
    }

    void MBEConfig::set_energy_threshold(double threshold) {
        if (threshold <= 0.0) {
            throw std::invalid_argument("Energy threshold must be positive");
        }
        energy_threshold_ = threshold;
    }

    void MBEConfig::set_qm_method(const std::string& method) {
        if (method.empty()) {
            throw std::invalid_argument("QM method cannot be empty");
        }
        qm_method_ = method;
    }

    void MBEConfig::set_basis_set(const std::string& basis) {
        if (basis.empty()) {
            throw std::invalid_argument("Basis set cannot be empty");
        }
        basis_set_ = basis;
    }

    void MBEConfig::set_charge_embedding(bool enable) {
        charge_embedding_ = enable;
    }

    void MBEConfig::set_custom_fragments(const std::vector<std::vector<std::size_t>>& fragments) {
        custom_fragments_ = fragments;
        fragmentation_scheme_ = FragmentationScheme::CUSTOM;
    }

    void MBEConfig::set_qm_options(const std::unordered_map<std::string, std::string>& options) {
        qm_options_ = options;
    }

    bool MBEConfig::is_valid() const {
        return validate_order() && validate_cutoffs() && validate_qm_settings();
    }

    std::string MBEConfig::to_string() const {
        std::ostringstream oss;
        oss << "MBEConfig:\n";
        oss << "  Max Order: " << max_order_ << "\n";
        oss << "  Fragmentation: " << libfrag::to_string(fragmentation_scheme_) << "\n";
        oss << "  Distance Cutoff: " << distance_cutoff_ << " Å\n";
        oss << "  Energy Threshold: " << energy_threshold_ << " Hartree\n";
        oss << "  QM Method: " << qm_method_ << "\n";
        oss << "  Basis Set: " << basis_set_ << "\n";
        oss << "  Charge Embedding: " << (charge_embedding_ ? "enabled" : "disabled") << "\n";
        
        if (!custom_fragments_.empty()) {
            oss << "  Custom Fragments: " << custom_fragments_.size() << " definitions\n";
        }
        
        if (!qm_options_.empty()) {
            oss << "  QM Options: " << qm_options_.size() << " options\n";
        }
        
        return oss.str();
    }

    MBEConfig MBEConfig::from_string(const std::string& config_string) {
        // TODO: Implement JSON/YAML parsing for configuration
        // This is a placeholder implementation
        throw std::runtime_error("Configuration parsing not yet implemented");
    }

    MBEConfig MBEConfig::default_2body() {
        MBEConfig config(2, FragmentationScheme::MOLECULAR);
        config.set_qm_method("HF");
        config.set_basis_set("6-31G*");
        config.set_distance_cutoff(8.0);
        config.set_energy_threshold(1e-6);
        return config;
    }

    MBEConfig MBEConfig::default_3body() {
        MBEConfig config(3, FragmentationScheme::MOLECULAR);
        config.set_qm_method("HF");
        config.set_basis_set("6-31G*");
        config.set_distance_cutoff(6.0);
        config.set_energy_threshold(1e-5);
        return config;
    }

    bool MBEConfig::validate_order() const {
        return max_order_ > 0 && max_order_ <= 10;  // Reasonable upper limit
    }

    bool MBEConfig::validate_cutoffs() const {
        return distance_cutoff_ > 0.0 && distance_cutoff_ < 100.0 &&
               energy_threshold_ > 0.0 && energy_threshold_ < 1.0;
    }

    bool MBEConfig::validate_qm_settings() const {
        return !qm_method_.empty() && !basis_set_.empty();
    }

    std::string to_string(MBEConfig::FragmentationScheme scheme) {
        switch (scheme) {
            case MBEConfig::FragmentationScheme::ATOMIC:
                return "atomic";
            case MBEConfig::FragmentationScheme::MOLECULAR:
                return "molecular";
            case MBEConfig::FragmentationScheme::FUNCTIONAL_GROUP:
                return "functional_group";
            case MBEConfig::FragmentationScheme::CUSTOM:
                return "custom";
            case MBEConfig::FragmentationScheme::DISTANCE_BASED:
                return "distance_based";
            default:
                return "unknown";
        }
    }

    std::string to_string(MBEConfig::TruncationMethod method) {
        switch (method) {
            case MBEConfig::TruncationMethod::ORDER_BASED:
                return "order_based";
            case MBEConfig::TruncationMethod::DISTANCE_BASED:
                return "distance_based";
            case MBEConfig::TruncationMethod::ENERGY_BASED:
                return "energy_based";
            default:
                return "unknown";
        }
    }

} // namespace libfrag
