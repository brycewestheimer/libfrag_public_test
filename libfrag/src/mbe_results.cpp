#include "libfrag/mbe_results.hpp"
#include "libfrag/mbe_config.hpp"
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <sstream>

namespace libfrag {

    // FragmentCalculationResult implementation
    FragmentCalculationResult::FragmentCalculationResult(
        const std::vector<std::size_t>& indices, int order)
        : fragment_indices(indices), n_body_order(order) {
        
        // Generate a default fragment ID
        std::ostringstream oss;
        oss << "frag_" << order << "body_";
        for (size_t i = 0; i < indices.size(); ++i) {
            if (i > 0) oss << "_";
            oss << indices[i];
        }
        fragment_id = oss.str();
    }

    std::string FragmentCalculationResult::to_string() const {
        std::ostringstream oss;
        oss << "Fragment " << fragment_id << " (" << n_body_order << "-body):\n";
        oss << "  Total Energy: " << std::fixed << std::setprecision(8) 
            << total_energy << " Hartree\n";
        oss << "  Method: " << qm_method << "/" << basis_set << "\n";
        oss << "  Converged: " << (converged ? "Yes" : "No") << "\n";
        oss << "  Time: " << computation_time.count() << " seconds\n";
        
        if (!properties.empty()) {
            oss << "  Properties:\n";
            for (const auto& [name, value] : properties) {
                oss << "    " << name << ": " << value << "\n";
            }
        }
        
        return oss.str();
    }

    bool FragmentCalculationResult::is_valid() const {
        return !fragment_indices.empty() && 
               n_body_order > 0 && 
               converged && 
               !fragment_id.empty();
    }

    // MBEResults implementation
    MBEResults::MBEResults(const MBEConfig& config) 
        : config_(std::make_shared<MBEConfig>(config)) {
        start_time_ = std::chrono::high_resolution_clock::now();
    }

    void MBEResults::add_fragment_result(const FragmentCalculationResult& result) {
        fragment_results_.push_back(result);
        max_order_ = std::max(max_order_, result.n_body_order);
        update_totals();
    }

    void MBEResults::add_fragment_result(FragmentCalculationResult&& result) {
        max_order_ = std::max(max_order_, result.n_body_order);
        fragment_results_.push_back(std::move(result));
        update_totals();
    }

    double MBEResults::energy_contribution(int order) const {
        auto it = energy_by_order_.find(order);
        return (it != energy_by_order_.end()) ? it->second : 0.0;
    }

    std::unordered_map<int, double> MBEResults::energy_contributions() const {
        return energy_by_order_;
    }

    std::vector<double> MBEResults::fragment_energies() const {
        std::vector<double> energies;
        for (const auto& result : fragment_results_) {
            if (result.n_body_order == 1) {
                energies.push_back(result.total_energy);
            }
        }
        return energies;
    }

    std::vector<double> MBEResults::interaction_energies(int order) const {
        std::vector<double> energies;
        for (const auto& result : fragment_results_) {
            if (result.n_body_order == order && order > 1) {
                // Calculate interaction energy (implementation placeholder)
                // This would involve subtracting lower-order contributions
                energies.push_back(result.total_energy);
            }
        }
        return energies;
    }

    std::vector<FragmentCalculationResult> MBEResults::results_by_order(int order) const {
        std::vector<FragmentCalculationResult> results;
        std::copy_if(fragment_results_.begin(), fragment_results_.end(),
                     std::back_inserter(results),
                     [order](const FragmentCalculationResult& r) {
                         return r.n_body_order == order;
                     });
        return results;
    }

    double MBEResults::estimated_truncation_error() const {
        if (max_order_ < 2) return 0.0;
        
        // Simple estimate based on energy contributions
        // More sophisticated methods would be implemented here
        double last_contribution = energy_contribution(max_order_);
        return std::abs(last_contribution) * 0.1;  // 10% of last term as rough estimate
    }

    bool MBEResults::is_converged(double threshold) const {
        if (max_order_ < 2) return false;
        
        double error = estimated_truncation_error();
        return error < threshold;
    }

    std::unordered_map<std::string, double> MBEResults::convergence_analysis() const {
        std::unordered_map<std::string, double> analysis;
        
        analysis["max_order"] = static_cast<double>(max_order_);
        analysis["total_energy"] = total_energy_;
        analysis["truncation_error"] = estimated_truncation_error();
        analysis["n_calculations"] = static_cast<double>(fragment_results_.size());
        
        // Energy contributions by order
        for (const auto& [order, energy] : energy_by_order_) {
            analysis["energy_" + std::to_string(order) + "body"] = energy;
        }
        
        return analysis;
    }

    std::chrono::duration<double> MBEResults::total_computation_time() const {
        if (end_time_.time_since_epoch().count() == 0) {
            // Calculation still ongoing
            auto now = std::chrono::high_resolution_clock::now();
            return now - start_time_;
        } else {
            return end_time_ - start_time_;
        }
    }

    std::chrono::duration<double> MBEResults::computation_time_by_order(int order) const {
        std::chrono::duration<double> total_time(0);
        for (const auto& result : fragment_results_) {
            if (result.n_body_order == order) {
                total_time += result.computation_time;
            }
        }
        return total_time;
    }

    std::unordered_map<std::string, double> MBEResults::performance_statistics() const {
        std::unordered_map<std::string, double> stats;
        
        auto total_time = total_computation_time();
        stats["total_time_seconds"] = total_time.count();
        stats["n_calculations"] = static_cast<double>(fragment_results_.size());
        
        if (!fragment_results_.empty()) {
            stats["avg_time_per_calc"] = total_time.count() / fragment_results_.size();
        }
        
        // Time by order
        for (int order = 1; order <= max_order_; ++order) {
            auto order_time = computation_time_by_order(order);
            stats["time_" + std::to_string(order) + "body"] = order_time.count();
            
            auto order_results = results_by_order(order);
            if (!order_results.empty()) {
                stats["avg_time_" + std::to_string(order) + "body"] = 
                    order_time.count() / order_results.size();
            }
        }
        
        return stats;
    }

    std::string MBEResults::to_json() const {
        // TODO: Implement proper JSON serialization
        std::ostringstream oss;
        oss << "{\n";
        oss << "  \"total_energy\": " << total_energy_ << ",\n";
        oss << "  \"max_order\": " << max_order_ << ",\n";
        oss << "  \"n_calculations\": " << fragment_results_.size() << ",\n";
        oss << "  \"energy_contributions\": {\n";
        
        bool first = true;
        for (const auto& [order, energy] : energy_by_order_) {
            if (!first) oss << ",\n";
            oss << "    \"" << order << "\": " << energy;
            first = false;
        }
        
        oss << "\n  }\n";
        oss << "}\n";
        return oss.str();
    }

    std::string MBEResults::to_csv() const {
        std::ostringstream oss;
        oss << "fragment_id,n_body_order,total_energy,converged,time_seconds\n";
        
        for (const auto& result : fragment_results_) {
            oss << result.fragment_id << ","
                << result.n_body_order << ","
                << std::fixed << std::setprecision(8) << result.total_energy << ","
                << (result.converged ? "true" : "false") << ","
                << result.computation_time.count() << "\n";
        }
        
        return oss.str();
    }

    std::string MBEResults::summary_report() const {
        std::ostringstream oss;
        oss << "=== MBE Calculation Summary ===\n\n";
        
        oss << "Total Energy: " << std::fixed << std::setprecision(8) 
            << total_energy_ << " Hartree\n";
        oss << "Maximum Order: " << max_order_ << "\n";
        oss << "Number of Calculations: " << fragment_results_.size() << "\n";
        oss << "Computation Time: " << total_computation_time().count() << " seconds\n\n";
        
        oss << "Energy Contributions by Order:\n";
        for (int order = 1; order <= max_order_; ++order) {
            double contribution = energy_contribution(order);
            oss << "  " << order << "-body: " << std::setw(12) << std::fixed 
                << std::setprecision(8) << contribution << " Hartree\n";
        }
        
        oss << "\nConvergence Analysis:\n";
        oss << "  Estimated Error: " << estimated_truncation_error() << " Hartree\n";
        oss << "  Converged (1e-6): " << (is_converged(1e-6) ? "Yes" : "No") << "\n";
        
        return oss.str();
    }

    std::string MBEResults::energy_decomposition_table() const {
        return format_energy_table();
    }

    bool MBEResults::validate() const {
        // Check basic consistency
        if (fragment_results_.empty()) return false;
        
        // Check that all results are valid
        for (const auto& result : fragment_results_) {
            if (!result.is_valid()) return false;
        }
        
        // Check energy consistency (placeholder)
        return true;
    }

    void MBEResults::clear() {
        fragment_results_.clear();
        total_energy_ = 0.0;
        max_order_ = 0;
        n_fragments_ = 0;
        energy_by_order_.clear();
    }

    const FragmentCalculationResult& MBEResults::operator[](std::size_t index) const {
        if (index >= fragment_results_.size()) {
            throw std::out_of_range("Index out of range");
        }
        return fragment_results_[index];
    }

    void MBEResults::update_totals() {
        compute_energy_contributions();
        
        // Update total energy (simple sum for now)
        total_energy_ = 0.0;
        for (const auto& [order, energy] : energy_by_order_) {
            total_energy_ += energy;
        }
        
        // Update fragment count estimate
        auto one_body_results = results_by_order(1);
        n_fragments_ = one_body_results.size();
    }

    void MBEResults::compute_energy_contributions() {
        energy_by_order_.clear();
        
        for (const auto& result : fragment_results_) {
            energy_by_order_[result.n_body_order] += result.total_energy;
        }
    }

    double MBEResults::calculate_interaction_energy(
        const std::vector<std::size_t>& fragment_indices) const {
        // TODO: Implement proper interaction energy calculation
        // This is a placeholder implementation
        return 0.0;
    }

    std::string MBEResults::format_energy_table() const {
        std::ostringstream oss;
        oss << std::setw(10) << "Order" << std::setw(15) << "Energy (Hartree)" 
            << std::setw(15) << "Contribution %" << "\n";
        oss << std::string(40, '-') << "\n";
        
        double total_abs = std::abs(total_energy_);
        
        for (int order = 1; order <= max_order_; ++order) {
            double contribution = energy_contribution(order);
            double percentage = (total_abs > 0) ? (contribution / total_abs) * 100.0 : 0.0;
            
            oss << std::setw(10) << order 
                << std::setw(15) << std::fixed << std::setprecision(8) << contribution
                << std::setw(14) << std::fixed << std::setprecision(2) << percentage << "%\n";
        }
        
        oss << std::string(40, '-') << "\n";
        oss << std::setw(10) << "Total" 
            << std::setw(15) << std::fixed << std::setprecision(8) << total_energy_ << "\n";
        
        return oss.str();
    }

    std::string MBEResults::format_timing_table() const {
        std::ostringstream oss;
        oss << std::setw(10) << "Order" << std::setw(15) << "Time (seconds)" 
            << std::setw(10) << "Count" << std::setw(15) << "Avg Time" << "\n";
        oss << std::string(50, '-') << "\n";
        
        for (int order = 1; order <= max_order_; ++order) {
            auto order_time = computation_time_by_order(order);
            auto order_results = results_by_order(order);
            double avg_time = order_results.empty() ? 0.0 : 
                              order_time.count() / order_results.size();
            
            oss << std::setw(10) << order 
                << std::setw(15) << std::fixed << std::setprecision(3) << order_time.count()
                << std::setw(10) << order_results.size()
                << std::setw(15) << std::fixed << std::setprecision(3) << avg_time << "\n";
        }
        
        return oss.str();
    }

    // Operators
    bool operator==(const FragmentCalculationResult& lhs, const FragmentCalculationResult& rhs) {
        return lhs.fragment_id == rhs.fragment_id &&
               lhs.n_body_order == rhs.n_body_order &&
               std::abs(lhs.total_energy - rhs.total_energy) < 1e-10;
    }

    std::ostream& operator<<(std::ostream& os, const FragmentCalculationResult& result) {
        os << result.to_string();
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const MBEResults& results) {
        os << results.summary_report();
        return os;
    }

} // namespace libfrag
