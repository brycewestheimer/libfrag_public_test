#include "libfrag/global_scf_results.hpp"
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iomanip>

namespace libfrag {

    // FragmentWaveFunction implementation
    FragmentWaveFunction::FragmentWaveFunction(std::size_t id, const std::vector<std::size_t>& atoms)
        : fragment_id(id), atom_indices(atoms) {}

    std::string FragmentWaveFunction::to_string() const {
        std::ostringstream oss;
        oss << "FragmentWaveFunction " << fragment_id << ":\n";
        oss << "  Name: " << fragment_name << "\n";
        oss << "  Atoms: " << atom_indices.size() << "\n";
        oss << "  Total energy: " << total_energy << " Eh\n";
        oss << "  Electronic energy: " << electronic_energy << " Eh\n";
        oss << "  Embedding energy: " << embedding_energy << " Eh\n";
        oss << "  Method: " << qm_method << "/" << basis_set << "\n";
        oss << "  Converged: " << (converged ? "yes" : "no") << "\n";
        oss << "  SCF iterations: " << scf_iterations << "\n";
        return oss.str();
    }

    bool FragmentWaveFunction::is_valid() const {
        return fragment_id != SIZE_MAX && 
               !atom_indices.empty() && 
               n_electrons > 0 && 
               multiplicity > 0;
    }

    double FragmentWaveFunction::get_property(const std::string& name, double default_value) const {
        auto it = properties.find(name);
        return (it != properties.end()) ? it->second : default_value;
    }

    void FragmentWaveFunction::set_property(const std::string& name, double value) {
        properties[name] = value;
    }

    double FragmentWaveFunction::electronic_charge() const {
        return -static_cast<double>(n_electrons);
    }

    double FragmentWaveFunction::nuclear_charge() const {
        // TODO: Calculate nuclear charge from atom types
        return static_cast<double>(n_electrons) + static_cast<double>(charge);
    }

    std::array<double, 3> FragmentWaveFunction::electric_dipole_moment() const {
        // TODO: Calculate dipole moment from density matrix and geometry
        return {0.0, 0.0, 0.0};
    }

    double FragmentWaveFunction::total_electron_density() const {
        // TODO: Calculate total electron density from density matrix
        return static_cast<double>(n_electrons);
    }

    // GlobalSCFIteration implementation
    GlobalSCFIteration::GlobalSCFIteration(int iter_num) : iteration_number(iter_num) {}

    std::string GlobalSCFIteration::to_string() const {
        std::ostringstream oss;
        oss << "SCF Iteration " << iteration_number << ":\n";
        oss << "  Total energy: " << std::fixed << std::setprecision(8) << total_energy << " Eh\n";
        oss << "  Energy change: " << std::scientific << std::setprecision(2) << energy_change << " Eh\n";
        oss << "  Max density change: " << std::scientific << std::setprecision(2) << max_density_change << "\n";
        oss << "  RMS density change: " << std::scientific << std::setprecision(2) << rms_density_change << "\n";
        oss << "  Converged: " << (overall_converged ? "yes" : "no") << "\n";
        if (diis_used) {
            oss << "  DIIS: dimension=" << diis_dimension << ", error=" << std::setprecision(2) << diis_error << "\n";
        }
        oss << "  Time: " << iteration_time.count() << " s\n";
        return oss.str();
    }

    bool GlobalSCFIteration::is_converged(const GlobalSCFConfig& config) const {
        bool energy_ok = std::abs(energy_change) < config.energy_threshold();
        bool density_ok = max_density_change < config.density_threshold();
        bool orbital_ok = max_orbital_change < config.orbital_threshold();
        bool polarization_ok = !config.polarization_enabled() || 
                              (polarization_change < config.polarization_threshold());

        switch (config.convergence_type()) {
            case GlobalSCFConfig::ConvergenceType::ENERGY:
                return energy_ok;
            case GlobalSCFConfig::ConvergenceType::DENSITY:
                return density_ok;
            case GlobalSCFConfig::ConvergenceType::ORBITAL:
                return orbital_ok;
            case GlobalSCFConfig::ConvergenceType::COMBINED:
                return energy_ok && density_ok && orbital_ok && polarization_ok;
            default:
                return false;
        }
    }

    void GlobalSCFIteration::set_diagnostic(const std::string& name, double value) {
        diagnostics[name] = value;
    }

    double GlobalSCFIteration::get_diagnostic(const std::string& name, double default_value) const {
        auto it = diagnostics.find(name);
        return (it != diagnostics.end()) ? it->second : default_value;
    }

    // GlobalSCFResults implementation
    GlobalSCFResults::GlobalSCFResults(const GlobalSCFConfig& config) : config_(config) {
        start_time_ = std::chrono::high_resolution_clock::now();
    }

    void GlobalSCFResults::add_fragment_wavefunction(FragmentWaveFunction wf) {
        std::size_t index = fragment_wavefunctions_.size();
        update_fragment_index(wf.fragment_id, index);
        fragment_wavefunctions_.push_back(std::move(wf));
    }

    const FragmentWaveFunction& GlobalSCFResults::get_fragment_wavefunction(std::size_t fragment_id) const {
        std::size_t index = get_fragment_index(fragment_id);
        if (index >= fragment_wavefunctions_.size()) {
            throw std::out_of_range("Fragment ID not found: " + std::to_string(fragment_id));
        }
        return fragment_wavefunctions_[index];
    }

    FragmentWaveFunction& GlobalSCFResults::get_fragment_wavefunction(std::size_t fragment_id) {
        std::size_t index = get_fragment_index(fragment_id);
        if (index >= fragment_wavefunctions_.size()) {
            throw std::out_of_range("Fragment ID not found: " + std::to_string(fragment_id));
        }
        return fragment_wavefunctions_[index];
    }

    void GlobalSCFResults::add_scf_iteration(GlobalSCFIteration iteration) {
        scf_iterations_.push_back(std::move(iteration));
    }

    const GlobalSCFIteration& GlobalSCFResults::get_scf_iteration(int iteration_number) const {
        if (iteration_number <= 0 || iteration_number > static_cast<int>(scf_iterations_.size())) {
            throw std::out_of_range("Invalid iteration number: " + std::to_string(iteration_number));
        }
        return scf_iterations_[iteration_number - 1];  // Convert to 0-based index
    }

    double GlobalSCFResults::final_energy() const {
        if (scf_iterations_.empty()) {
            return 0.0;
        }
        return scf_iterations_.back().total_energy;
    }

    std::unordered_map<std::string, double> GlobalSCFResults::energy_breakdown() const {
        std::unordered_map<std::string, double> breakdown;
        if (scf_iterations_.empty()) {
            return breakdown;
        }

        const auto& final_iter = scf_iterations_.back();
        breakdown["total_energy"] = final_iter.total_energy;
        breakdown["electronic_energy"] = final_iter.electronic_energy;
        breakdown["nuclear_repulsion"] = final_iter.nuclear_repulsion;
        breakdown["embedding_energy"] = final_iter.embedding_energy;
        breakdown["polarization_energy"] = final_iter.polarization_energy;
        breakdown["interaction_energy"] = final_iter.interaction_energy;

        return breakdown;
    }

    std::unordered_map<std::size_t, double> GlobalSCFResults::fragment_energies() const {
        std::unordered_map<std::size_t, double> energies;
        for (const auto& wf : fragment_wavefunctions_) {
            energies[wf.fragment_id] = wf.total_energy;
        }
        return energies;
    }

    double GlobalSCFResults::interaction_energy() const {
        if (scf_iterations_.empty()) {
            return 0.0;
        }
        return scf_iterations_.back().interaction_energy;
    }

    bool GlobalSCFResults::is_converged() const {
        if (scf_iterations_.empty()) {
            return false;
        }
        return scf_iterations_.back().overall_converged;
    }

    std::unordered_map<std::string, double> GlobalSCFResults::convergence_analysis() const {
        std::unordered_map<std::string, double> analysis;
        
        if (scf_iterations_.empty()) {
            return analysis;
        }

        const auto& final_iter = scf_iterations_.back();
        analysis["final_energy_change"] = final_iter.energy_change;
        analysis["final_density_change"] = final_iter.max_density_change;
        analysis["final_orbital_change"] = final_iter.max_orbital_change;
        analysis["iterations_to_convergence"] = static_cast<double>(scf_iterations_.size());
        analysis["converged"] = is_converged() ? 1.0 : 0.0;

        if (scf_iterations_.size() > 1) {
            // Calculate convergence rate
            double initial_energy_change = scf_iterations_[0].energy_change;
            double final_energy_change = final_iter.energy_change;
            if (initial_energy_change != 0.0) {
                analysis["convergence_rate"] = std::log10(std::abs(final_energy_change / initial_energy_change)) / 
                                              static_cast<double>(scf_iterations_.size() - 1);
            }
        }

        return analysis;
    }

    std::vector<double> GlobalSCFResults::convergence_history(const std::string& criterion) const {
        std::vector<double> history;
        history.reserve(scf_iterations_.size());

        for (const auto& iteration : scf_iterations_) {
            if (criterion == "energy_change") {
                history.push_back(std::abs(iteration.energy_change));
            } else if (criterion == "density_change") {
                history.push_back(iteration.max_density_change);
            } else if (criterion == "orbital_change") {
                history.push_back(iteration.max_orbital_change);
            } else if (criterion == "polarization_change") {
                history.push_back(iteration.polarization_change);
            } else if (criterion == "total_energy") {
                history.push_back(iteration.total_energy);
            } else {
                // Try diagnostics
                history.push_back(iteration.get_diagnostic(criterion, 0.0));
            }
        }

        return history;
    }

    std::chrono::duration<double> GlobalSCFResults::total_time() const {
        if (scf_iterations_.empty()) {
            return std::chrono::duration<double>(0.0);
        }

        return std::accumulate(scf_iterations_.begin(), scf_iterations_.end(),
                              std::chrono::duration<double>(0.0),
                              [](const std::chrono::duration<double>& sum, const GlobalSCFIteration& iter) {
                                  return sum + iter.iteration_time;
                              });
    }

    std::unordered_map<std::string, double> GlobalSCFResults::performance_breakdown() const {
        std::unordered_map<std::string, double> breakdown;
        
        double total_iteration_time = 0.0;
        double total_fragment_time = 0.0;
        double total_embedding_time = 0.0;

        for (const auto& iteration : scf_iterations_) {
            total_iteration_time += iteration.iteration_time.count();
            total_fragment_time += iteration.fragment_time.count();
            total_embedding_time += iteration.embedding_time.count();
        }

        breakdown["total_time"] = total_iteration_time;
        breakdown["fragment_time"] = total_fragment_time;
        breakdown["embedding_time"] = total_embedding_time;
        breakdown["overhead_time"] = total_iteration_time - total_fragment_time - total_embedding_time;

        if (total_iteration_time > 0.0) {
            breakdown["fragment_fraction"] = total_fragment_time / total_iteration_time;
            breakdown["embedding_fraction"] = total_embedding_time / total_iteration_time;
            breakdown["overhead_fraction"] = breakdown["overhead_time"] / total_iteration_time;
        }

        return breakdown;
    }

    std::unordered_map<std::string, double> GlobalSCFResults::performance_statistics() const {
        std::unordered_map<std::string, double> stats;
        
        if (scf_iterations_.empty()) {
            return stats;
        }

        std::vector<double> iteration_times;
        for (const auto& iteration : scf_iterations_) {
            iteration_times.push_back(iteration.iteration_time.count());
        }

        stats["avg_iteration_time"] = std::accumulate(iteration_times.begin(), iteration_times.end(), 0.0) / 
                                     iteration_times.size();
        
        auto minmax = std::minmax_element(iteration_times.begin(), iteration_times.end());
        stats["min_iteration_time"] = *minmax.first;
        stats["max_iteration_time"] = *minmax.second;

        stats["total_iterations"] = static_cast<double>(scf_iterations_.size());
        stats["time_per_fragment"] = stats["avg_iteration_time"] / static_cast<double>(fragment_wavefunctions_.size());

        return stats;
    }

    void GlobalSCFResults::export_json(const std::string& filename, bool include_wavefunctions) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file for writing: " + filename);
        }
        file << export_json_string(include_wavefunctions);
    }

    void GlobalSCFResults::export_convergence_csv(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file for writing: " + filename);
        }

        // Write header
        std::vector<std::string> headers = {
            "iteration", "total_energy", "energy_change", "max_density_change", 
            "rms_density_change", "max_orbital_change", "polarization_change",
            "converged", "iteration_time"
        };
        write_csv_header(file, headers);

        // Write data
        for (const auto& iteration : scf_iterations_) {
            std::vector<double> values = {
                static_cast<double>(iteration.iteration_number),
                iteration.total_energy,
                iteration.energy_change,
                iteration.max_density_change,
                iteration.rms_density_change,
                iteration.max_orbital_change,
                iteration.polarization_change,
                iteration.overall_converged ? 1.0 : 0.0,
                iteration.iteration_time.count()
            };
            write_csv_row(file, values);
        }
    }

    void GlobalSCFResults::export_energy_csv(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file for writing: " + filename);
        }

        // Write header
        std::vector<std::string> headers = {
            "fragment_id", "fragment_name", "total_energy", "electronic_energy",
            "nuclear_repulsion", "embedding_energy", "polarization_energy"
        };
        write_csv_header(file, headers);

        // Write fragment data
        for (const auto& wf : fragment_wavefunctions_) {
            file << wf.fragment_id << ","
                 << "\"" << wf.fragment_name << "\","
                 << wf.total_energy << ","
                 << wf.electronic_energy << ","
                 << wf.nuclear_repulsion << ","
                 << wf.embedding_energy << ","
                 << wf.polarization_energy << "\n";
        }
    }

    std::string GlobalSCFResults::generate_summary(bool detailed) const {
        std::ostringstream oss;
        
        oss << "=== Global SCF Results Summary ===\n\n";
        
        // Basic information
        oss << "Configuration:\n";
        oss << "  Method: " << config_.qm_method() << "/" << config_.basis_set() << "\n";
        oss << "  Embedding: " << GlobalSCFConfig::embedding_type_to_string(config_.embedding_type()) << "\n";
        oss << "  Fragmentation: " << GlobalSCFConfig::fragmentation_scheme_to_string(config_.fragmentation_scheme()) << "\n";
        oss << "\n";

        // Results summary
        oss << "Results:\n";
        oss << "  Final energy: " << std::fixed << std::setprecision(8) << final_energy() << " Eh\n";
        oss << "  Converged: " << (is_converged() ? "yes" : "no") << "\n";
        oss << "  SCF iterations: " << scf_iterations_.size() << "\n";
        oss << "  Number of fragments: " << fragment_wavefunctions_.size() << "\n";
        oss << "  Total time: " << total_time().count() << " s\n";

        if (detailed && !scf_iterations_.empty()) {
            oss << "\nConvergence Analysis:\n";
            auto conv_analysis = convergence_analysis();
            for (const auto& [key, value] : conv_analysis) {
                oss << "  " << key << ": " << std::scientific << std::setprecision(2) << value << "\n";
            }

            oss << "\nEnergy Breakdown:\n";
            auto energy_bd = energy_breakdown();
            for (const auto& [component, energy] : energy_bd) {
                oss << "  " << component << ": " << std::fixed << std::setprecision(6) << energy << " Eh\n";
            }

            oss << "\nFragment Energies:\n";
            auto frag_energies = fragment_energies();
            for (const auto& [frag_id, energy] : frag_energies) {
                oss << "  Fragment " << frag_id << ": " << std::fixed << std::setprecision(6) << energy << " Eh\n";
            }
        }

        return oss.str();
    }

    std::vector<std::string> GlobalSCFResults::validate() const {
        std::vector<std::string> issues;

        // Check basic consistency
        if (fragment_wavefunctions_.empty()) {
            issues.push_back("No fragment wave functions stored");
        }

        if (scf_iterations_.empty()) {
            issues.push_back("No SCF iterations recorded");
        }

        // Check fragment wave function validity
        for (const auto& wf : fragment_wavefunctions_) {
            if (!wf.is_valid()) {
                issues.push_back("Invalid fragment wave function: " + std::to_string(wf.fragment_id));
            }
        }

        // Check iteration consistency
        for (std::size_t i = 0; i < scf_iterations_.size(); ++i) {
            const auto& iter = scf_iterations_[i];
            if (iter.iteration_number != static_cast<int>(i) + 1) {
                issues.push_back("Inconsistent iteration numbering at iteration " + std::to_string(i + 1));
            }
        }

        return issues;
    }

    std::unordered_map<std::string, double> GlobalSCFResults::calculation_statistics() const {
        std::unordered_map<std::string, double> stats;
        
        stats["n_fragments"] = static_cast<double>(fragment_wavefunctions_.size());
        stats["n_iterations"] = static_cast<double>(scf_iterations_.size());
        stats["final_energy"] = final_energy();
        stats["converged"] = is_converged() ? 1.0 : 0.0;
        stats["total_time"] = total_time().count();

        if (!fragment_wavefunctions_.empty()) {
            // Fragment statistics
            std::vector<std::size_t> fragment_sizes;
            int total_electrons = 0;
            for (const auto& wf : fragment_wavefunctions_) {
                fragment_sizes.push_back(wf.atom_indices.size());
                total_electrons += wf.n_electrons;
            }

            auto minmax_size = std::minmax_element(fragment_sizes.begin(), fragment_sizes.end());
            stats["min_fragment_size"] = static_cast<double>(*minmax_size.first);
            stats["max_fragment_size"] = static_cast<double>(*minmax_size.second);
            stats["avg_fragment_size"] = static_cast<double>(std::accumulate(fragment_sizes.begin(), fragment_sizes.end(), 0.0)) / 
                                       static_cast<double>(fragment_sizes.size());
            stats["total_electrons"] = static_cast<double>(total_electrons);
        }

        return stats;
    }

    void GlobalSCFResults::clear() {
        fragment_wavefunctions_.clear();
        fragment_id_to_index_.clear();
        scf_iterations_.clear();
        start_time_ = std::chrono::high_resolution_clock::now();
        end_time_ = start_time_;
    }

    std::string GlobalSCFResults::to_string() const {
        return generate_summary(false);
    }

    // Private utility methods
    void GlobalSCFResults::update_fragment_index(std::size_t fragment_id, std::size_t index) {
        fragment_id_to_index_[fragment_id] = index;
    }

    std::size_t GlobalSCFResults::get_fragment_index(std::size_t fragment_id) const {
        auto it = fragment_id_to_index_.find(fragment_id);
        if (it == fragment_id_to_index_.end()) {
            return SIZE_MAX;
        }
        return it->second;
    }

    std::string GlobalSCFResults::export_json_string(bool include_wavefunctions) const {
        // TODO: Implement JSON export
        // For now, return a placeholder
        return "{\n  \"message\": \"JSON export not yet implemented\"\n}";
    }

    void GlobalSCFResults::write_csv_header(std::ostream& os, const std::vector<std::string>& headers) const {
        for (std::size_t i = 0; i < headers.size(); ++i) {
            if (i > 0) os << ",";
            os << headers[i];
        }
        os << "\n";
    }

    void GlobalSCFResults::write_csv_row(std::ostream& os, const std::vector<double>& values) const {
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i > 0) os << ",";
            os << std::scientific << std::setprecision(6) << values[i];
        }
        os << "\n";
    }

}
