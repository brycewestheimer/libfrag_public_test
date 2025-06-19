#include "libfrag/mbe_calculator.hpp"
#include <thread>
#include <algorithm>
#include <future>
#include <stdexcept>

namespace libfrag {

    // MBECalculator implementation
    MBECalculator::MBECalculator(QMCalculatorPtr qm_calculator) 
        : qm_calculator_(std::move(qm_calculator)) {}

    MBECalculator::MBECalculator(QMCalculatorPtr qm_calculator, const MBEConfig& config)
        : config_(config), qm_calculator_(std::move(qm_calculator)) {}

    void MBECalculator::set_qm_calculator(QMCalculatorPtr qm_calculator) {
        qm_calculator_ = std::move(qm_calculator);
    }

    MBEResults MBECalculator::calculate(const Molecule& molecule) {
        return calculate(molecule, config_);
    }

    MBEResults MBECalculator::calculate(const Molecule& molecule, const MBEConfig& config) {
        config_ = config;
        
        // Validate setup
        auto validation_errors = validate_setup(molecule, config);
        if (!validation_errors.empty()) {
            std::string error_msg = "Validation failed:\n";
            for (const auto& error : validation_errors) {
                error_msg += "  " + error + "\n";
            }
            throw std::runtime_error(error_msg);
        }
        
        update_progress(0.0, "Starting MBE calculation");
        
        return perform_calculation(molecule);
    }

    std::future<MBEResults> MBECalculator::calculate_async(const Molecule& molecule) {
        return std::async(std::launch::async, [this, molecule]() {
            return calculate(molecule);
        });
    }

    std::string MBECalculator::start_calculation(const Molecule& molecule) {
        // Generate unique calculation ID
        std::string calc_id = "mbe_calc_" + std::to_string(
            std::chrono::high_resolution_clock::now().time_since_epoch().count());
        
        // Start async calculation
        auto future = std::async(std::launch::async, [this, molecule]() {
            return calculate(molecule);
        });
        
        ongoing_calculations_[calc_id] = std::move(future);
        
        return calc_id;
    }

    MBEResults MBECalculator::get_partial_results(const std::string& calculation_id) {
        // TODO: Implement partial results tracking
        // For now, return empty results
        return MBEResults();
    }

    bool MBECalculator::is_calculation_complete(const std::string& calculation_id) {
        auto it = ongoing_calculations_.find(calculation_id);
        if (it == ongoing_calculations_.end()) {
            return false;
        }
        
        return it->second.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    }

    void MBECalculator::cancel_calculation(const std::string& calculation_id) {
        // TODO: Implement calculation cancellation
        // For now, just remove from tracking
        ongoing_calculations_.erase(calculation_id);
    }

    std::unordered_map<std::string, double> MBECalculator::estimate_cost(
        const Molecule& molecule, const MBEConfig& config) {
        
        std::unordered_map<std::string, double> cost_estimate;
        
        // Generate fragment combinations to estimate cost
        fragment_generator_.set_config(config);
        auto all_combinations = fragment_generator_.generate_all_combinations(molecule, config);
        
        std::size_t total_calculations = 0;
        for (const auto& [order, combinations] : all_combinations) {
            total_calculations += combinations.size();
        }
        
        cost_estimate["total_calculations"] = static_cast<double>(total_calculations);
        
        // Rough time estimates (these would be calibrated based on actual performance)
        double estimated_seconds = total_calculations * 10.0;  // 10 seconds per calculation
        cost_estimate["estimated_time_seconds"] = estimated_seconds;
        cost_estimate["estimated_time_hours"] = estimated_seconds / 3600.0;
        
        // Memory estimates
        double estimated_memory_mb = total_calculations * 100.0;  // 100 MB per calculation
        cost_estimate["estimated_memory_mb"] = estimated_memory_mb;
        
        return cost_estimate;
    }

    std::unordered_map<int, std::size_t> MBECalculator::preview_fragments(
        const Molecule& molecule, const MBEConfig& config) {
        
        fragment_generator_.set_config(config);
        auto all_combinations = fragment_generator_.generate_all_combinations(molecule, config);
        
        std::unordered_map<int, std::size_t> fragment_counts;
        for (const auto& [order, combinations] : all_combinations) {
            fragment_counts[order] = combinations.size();
        }
        
        return fragment_counts;
    }

    std::vector<std::string> MBECalculator::validate_setup(
        const Molecule& molecule, const MBEConfig& config) {
        
        std::vector<std::string> errors;
        
        // Validate molecule
        auto molecule_errors = validate_molecule(molecule);
        errors.insert(errors.end(), molecule_errors.begin(), molecule_errors.end());
        
        // Validate config
        auto config_errors = validate_config(config);
        errors.insert(errors.end(), config_errors.begin(), config_errors.end());
        
        // Validate QM calculator
        auto qm_errors = validate_qm_calculator();
        errors.insert(errors.end(), qm_errors.begin(), qm_errors.end());
        
        return errors;
    }

    bool MBECalculator::is_ready() const {
        return qm_calculator_ && qm_calculator_->is_available() && config_.is_valid();
    }

    std::unordered_map<std::string, std::string> MBECalculator::get_status() const {
        std::unordered_map<std::string, std::string> status;
        
        status["ready"] = is_ready() ? "true" : "false";
        status["qm_available"] = (qm_calculator_ && qm_calculator_->is_available()) ? "true" : "false";
        status["config_valid"] = config_.is_valid() ? "true" : "false";
        status["n_threads"] = std::to_string(n_threads_);
        status["caching_enabled"] = caching_enabled_ ? "true" : "false";
        status["ongoing_calculations"] = std::to_string(ongoing_calculations_.size());
        
        if (qm_calculator_) {
            status["qm_software"] = qm_calculator_->software_info();
        }
        
        return status;
    }

    void MBECalculator::clear_cache() {
        calculation_cache_.clear();
    }

    std::unordered_map<std::string, double> MBECalculator::get_performance_stats() const {
        return performance_stats_;
    }

    // Private methods
    MBEResults MBECalculator::perform_calculation(const Molecule& molecule) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        MBEResults results(config_);
        
        // Generate all fragment combinations
        fragment_generator_.set_config(config_);
        auto all_combinations = fragment_generator_.generate_all_combinations(molecule, config_);
        
        // Create initial fragments
        auto base_fragments = fragment_generator_.create_initial_fragments(molecule, config_);
        
        update_progress(10.0, "Generated fragment combinations");
        
        std::size_t total_calculations = 0;
        for (const auto& [order, combinations] : all_combinations) {
            total_calculations += combinations.size();
        }
        
        std::size_t completed_calculations = 0;
        
        // Perform calculations for each N-body order
        for (int order = 1; order <= config_.max_order(); ++order) {
            log_message("Calculating " + std::to_string(order) + "-body terms");
            
            auto combinations_it = all_combinations.find(order);
            if (combinations_it == all_combinations.end()) continue;
            
            const auto& combinations = combinations_it->second;
            
            // Calculate all combinations for this order
            for (const auto& combination : combinations) {
                try {
                    auto fragment_results = calculate_single_fragment(combination, base_fragments);
                    results.add_fragment_result(std::move(fragment_results));
                    
                    completed_calculations++;
                    double progress = 10.0 + (90.0 * completed_calculations / total_calculations);
                    update_progress(progress, "Completed " + std::to_string(completed_calculations) + 
                                   "/" + std::to_string(total_calculations) + " calculations");
                    
                } catch (const std::exception& e) {
                    handle_calculation_error(e, combination);
                }
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
        
        performance_stats_["last_calculation_time"] = duration.count();
        performance_stats_["last_calculation_count"] = static_cast<double>(total_calculations);
        
        update_progress(100.0, "MBE calculation complete");
        
        return results;
    }

    std::vector<FragmentCalculationResult> MBECalculator::calculate_n_body_terms(
        const std::vector<std::size_t>& fragment_indices,
        const std::vector<std::shared_ptr<Fragment>>& fragments,
        int n_body) {
        
        // TODO: Implement parallel calculation of N-body terms
        std::vector<FragmentCalculationResult> results;
        
        return results;
    }

    FragmentCalculationResult MBECalculator::calculate_single_fragment(
        const std::vector<std::size_t>& fragment_indices,
        const std::vector<std::shared_ptr<Fragment>>& fragments) {
        
        if (!qm_calculator_) {
            throw std::runtime_error("No QM calculator available");
        }
        
        // Create combined fragment from indices
        auto combined_fragment = fragment_generator_.create_combined_fragment(
            fragment_indices, fragments, 
            fragment_generator_.generate_combination_id(fragment_indices, fragment_indices.size()));
        
        if (!combined_fragment) {
            throw std::runtime_error("Failed to create combined fragment");
        }
        
        // Check cache first
        std::string cache_key = generate_cache_key(*combined_fragment, 
            config_.qm_method(), config_.basis_set(), 0, 1);
        
        if (caching_enabled_ && is_calculation_cached(cache_key)) {
            return get_cached_result(cache_key);
        }
        
        // Perform QM calculation
        auto result = qm_calculator_->calculate_energy(
            *combined_fragment,
            config_.qm_method(),
            config_.basis_set(),
            0,  // charge
            1,  // multiplicity
            config_.qm_options()
        );
        
        // Cache result
        if (caching_enabled_) {
            cache_result(cache_key, result);
        }
        
        return result;
    }

    void MBECalculator::update_progress(double progress, const std::string& status) {
        if (progress_callback_) {
            progress_callback_(progress, status);
        }
    }

    void MBECalculator::log_message(const std::string& message) {
        if (log_callback_) {
            log_callback_(message);
        }
    }

    std::string MBECalculator::generate_cache_key(const Fragment& fragment, 
        const std::string& method, const std::string& basis_set,
        int charge, int multiplicity) {
        
        // TODO: Generate a proper cache key based on fragment geometry and settings
        return fragment.fragment_id() + "_" + method + "_" + basis_set + 
               "_" + std::to_string(charge) + "_" + std::to_string(multiplicity);
    }

    bool MBECalculator::is_calculation_cached(const std::string& cache_key) {
        return calculation_cache_.find(cache_key) != calculation_cache_.end();
    }

    FragmentCalculationResult MBECalculator::get_cached_result(const std::string& cache_key) {
        auto it = calculation_cache_.find(cache_key);
        if (it != calculation_cache_.end()) {
            return it->second;
        }
        throw std::runtime_error("Cache key not found: " + cache_key);
    }

    void MBECalculator::cache_result(const std::string& cache_key, 
        const FragmentCalculationResult& result) {
        calculation_cache_[cache_key] = result;
    }

    std::vector<std::string> MBECalculator::validate_molecule(const Molecule& molecule) {
        std::vector<std::string> errors;
        
        if (molecule.size() == 0) {
            errors.push_back("Molecule is empty");
        }
        
        // TODO: Add more molecule validation
        
        return errors;
    }

    std::vector<std::string> MBECalculator::validate_config(const MBEConfig& config) {
        std::vector<std::string> errors;
        
        if (!config.is_valid()) {
            errors.push_back("Configuration is invalid");
        }
        
        return errors;
    }

    std::vector<std::string> MBECalculator::validate_qm_calculator() {
        std::vector<std::string> errors;
        
        if (!qm_calculator_) {
            errors.push_back("No QM calculator provided");
        } else if (!qm_calculator_->is_available()) {
            errors.push_back("QM calculator is not available");
        }
        
        return errors;
    }

    void MBECalculator::handle_calculation_error(const std::exception& e, 
        const std::vector<std::size_t>& fragment_indices) {
        
        std::string error_msg = "Calculation failed for fragment indices [";
        for (size_t i = 0; i < fragment_indices.size(); ++i) {
            if (i > 0) error_msg += ", ";
            error_msg += std::to_string(fragment_indices[i]);
        }
        error_msg += "]: " + std::string(e.what());
        
        log_message(error_msg);
        
        // TODO: Decide whether to continue or abort the calculation
    }

    // Factory and utility functions
    std::unique_ptr<MBECalculator> create_mbe_calculator(
        const std::string& qm_software,
        const std::unordered_map<std::string, std::string>& options) {
        
        // TODO: Implement QM software-specific calculator creation
        // This would create appropriate QMCalculatorInterface implementations
        
        throw std::runtime_error("QM software integration not yet implemented: " + qm_software);
    }

    MBEResults calculate_mbe(const Molecule& molecule, int max_order,
        const std::string& qm_method, const std::string& basis_set) {
        
        // TODO: Implement convenience function
        // This would create a default calculator and perform the calculation
        
        throw std::runtime_error("Convenience function not yet implemented");
    }

} // namespace libfrag
