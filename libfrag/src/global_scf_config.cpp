#include "libfrag/global_scf_config.hpp"
#include <sstream>
#include <stdexcept>

namespace libfrag {

    // GlobalSCFConfig implementation
    GlobalSCFConfig::GlobalSCFConfig(FragmentationScheme scheme, EmbeddingType embedding)
        : fragmentation_scheme_(scheme), embedding_type_(embedding) {}

    void GlobalSCFConfig::set_fragmentation_scheme(FragmentationScheme scheme) {
        fragmentation_scheme_ = scheme;
    }

    void GlobalSCFConfig::set_embedding_type(EmbeddingType embedding) {
        embedding_type_ = embedding;
    }

    void GlobalSCFConfig::set_convergence_type(ConvergenceType convergence) {
        convergence_type_ = convergence;
    }

    void GlobalSCFConfig::set_max_scf_iterations(int max_iter) {
        if (max_iter <= 0) {
            throw std::invalid_argument("Maximum SCF iterations must be positive");
        }
        max_scf_iterations_ = max_iter;
    }

    void GlobalSCFConfig::set_energy_threshold(double threshold) {
        if (threshold <= 0.0) {
            throw std::invalid_argument("Energy threshold must be positive");
        }
        energy_threshold_ = threshold;
    }

    void GlobalSCFConfig::set_density_threshold(double threshold) {
        if (threshold <= 0.0) {
            throw std::invalid_argument("Density threshold must be positive");
        }
        density_threshold_ = threshold;
    }

    void GlobalSCFConfig::set_orbital_threshold(double threshold) {
        if (threshold <= 0.0) {
            throw std::invalid_argument("Orbital threshold must be positive");
        }
        orbital_threshold_ = threshold;
    }

    void GlobalSCFConfig::set_qm_method(const std::string& method) {
        if (method.empty()) {
            throw std::invalid_argument("QM method cannot be empty");
        }
        qm_method_ = method;
    }

    void GlobalSCFConfig::set_basis_set(const std::string& basis) {
        if (basis.empty()) {
            throw std::invalid_argument("Basis set cannot be empty");
        }
        basis_set_ = basis;
    }

    void GlobalSCFConfig::set_diis_enabled(bool enable) {
        diis_enabled_ = enable;
    }

    void GlobalSCFConfig::set_diis_size(int size) {
        if (size < 2) {
            throw std::invalid_argument("DIIS size must be at least 2");
        }
        diis_size_ = size;
    }

    void GlobalSCFConfig::set_buffer_distance(double distance) {
        if (distance < 0.0) {
            throw std::invalid_argument("Buffer distance cannot be negative");
        }
        buffer_distance_ = distance;
    }

    void GlobalSCFConfig::set_custom_fragments(const std::vector<std::vector<std::size_t>>& fragments) {
        custom_fragments_ = fragments;
        fragmentation_scheme_ = FragmentationScheme::CUSTOM;
    }

    void GlobalSCFConfig::set_qm_options(const std::unordered_map<std::string, std::string>& options) {
        qm_options_ = options;
    }

    void GlobalSCFConfig::set_polarization_enabled(bool enable) {
        polarization_enabled_ = enable;
    }

    void GlobalSCFConfig::set_polarization_threshold(double threshold) {
        if (threshold <= 0.0) {
            throw std::invalid_argument("Polarization threshold must be positive");
        }
        polarization_threshold_ = threshold;
    }

    bool GlobalSCFConfig::is_valid() const {
        return validate_thresholds() && validate_qm_settings() && validate_diis_settings();
    }

    std::string GlobalSCFConfig::to_string() const {
        std::ostringstream oss;
        oss << "GlobalSCFConfig:\n";
        oss << "  Fragmentation: " << fragmentation_scheme_to_string(fragmentation_scheme_) << "\n";
        oss << "  Embedding: " << embedding_type_to_string(embedding_type_) << "\n";
        oss << "  Convergence: " << convergence_type_to_string(convergence_type_) << "\n";
        oss << "  Max SCF iterations: " << max_scf_iterations_ << "\n";
        oss << "  Energy threshold: " << energy_threshold_ << "\n";
        oss << "  Density threshold: " << density_threshold_ << "\n";
        oss << "  QM method: " << qm_method_ << "\n";
        oss << "  Basis set: " << basis_set_ << "\n";
        oss << "  DIIS enabled: " << (diis_enabled_ ? "yes" : "no") << "\n";
        oss << "  Buffer distance: " << buffer_distance_ << " Å\n";
        oss << "  Polarization enabled: " << (polarization_enabled_ ? "yes" : "no") << "\n";
        return oss.str();
    }

    // Factory methods
    GlobalSCFConfig GlobalSCFConfig::create_fmo_config(const std::string& method, const std::string& basis) {
        GlobalSCFConfig config(FragmentationScheme::FMO_LIKE, EmbeddingType::COULOMB);
        config.set_qm_method(method);
        config.set_basis_set(basis);
        config.set_diis_enabled(true);
        config.set_polarization_enabled(true);
        config.set_buffer_distance(2.0);
        return config;
    }

    GlobalSCFConfig GlobalSCFConfig::create_coulomb_config(const std::string& method, const std::string& basis) {
        GlobalSCFConfig config(FragmentationScheme::MOLECULAR, EmbeddingType::COULOMB);
        config.set_qm_method(method);
        config.set_basis_set(basis);
        config.set_diis_enabled(true);
        config.set_polarization_enabled(false);
        return config;
    }

    GlobalSCFConfig GlobalSCFConfig::create_polarizable_config(const std::string& method, const std::string& basis) {
        GlobalSCFConfig config(FragmentationScheme::MOLECULAR, EmbeddingType::POLARIZABLE);
        config.set_qm_method(method);
        config.set_basis_set(basis);
        config.set_diis_enabled(true);
        config.set_polarization_enabled(true);
        config.set_polarization_threshold(1e-6);
        return config;
    }

    // String conversion utilities
    std::string GlobalSCFConfig::fragmentation_scheme_to_string(FragmentationScheme scheme) {
        switch (scheme) {
            case FragmentationScheme::ATOMIC: return "atomic";
            case FragmentationScheme::MOLECULAR: return "molecular";
            case FragmentationScheme::FUNCTIONAL_GROUP: return "functional_group";
            case FragmentationScheme::CUSTOM: return "custom";
            case FragmentationScheme::DISTANCE_BASED: return "distance_based";
            case FragmentationScheme::FMO_LIKE: return "fmo_like";
            default: return "unknown";
        }
    }

    std::string GlobalSCFConfig::embedding_type_to_string(EmbeddingType embedding) {
        switch (embedding) {
            case EmbeddingType::NONE: return "none";
            case EmbeddingType::COULOMB: return "coulomb";
            case EmbeddingType::POLARIZABLE: return "polarizable";
            case EmbeddingType::DENSITY_BASED: return "density_based";
            case EmbeddingType::COMBINED: return "combined";
            default: return "unknown";
        }
    }

    std::string GlobalSCFConfig::convergence_type_to_string(ConvergenceType convergence) {
        switch (convergence) {
            case ConvergenceType::ENERGY: return "energy";
            case ConvergenceType::DENSITY: return "density";
            case ConvergenceType::ORBITAL: return "orbital";
            case ConvergenceType::COMBINED: return "combined";
            default: return "unknown";
        }
    }

    // Private validation methods
    bool GlobalSCFConfig::validate_thresholds() const {
        return energy_threshold_ > 0.0 && 
               density_threshold_ > 0.0 && 
               orbital_threshold_ > 0.0 &&
               polarization_threshold_ > 0.0;
    }

    bool GlobalSCFConfig::validate_qm_settings() const {
        return !qm_method_.empty() && !basis_set_.empty();
    }

    bool GlobalSCFConfig::validate_diis_settings() const {
        return !diis_enabled_ || diis_size_ >= 2;
    }

}
