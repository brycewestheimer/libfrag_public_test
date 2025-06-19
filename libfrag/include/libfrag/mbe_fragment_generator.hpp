#pragma once
#ifndef LIBFRAG_MBE_FRAGMENT_GENERATOR_HPP
#define LIBFRAG_MBE_FRAGMENT_GENERATOR_HPP

#include "libfrag/fragment.hpp"
#include "libfrag/molecule.hpp"
#include "libfrag/mbe_config.hpp"
#include <vector>
#include <memory>
#include <unordered_set>
#include <functional>

namespace libfrag {

    /**
     * @brief Generates fragment combinations for Many-Body Expansion
     * 
     * This class is responsible for creating all the necessary fragment combinations
     * for MBE calculations based on the specified fragmentation scheme and truncation
     * criteria. It handles the combinatorial generation of 1-body, 2-body, 3-body, etc.
     * fragment combinations.
     */
    class MBEFragmentGenerator {
    public:
        // Type aliases
        using FragmentCombination = std::vector<std::size_t>;
        using FragmentCombinations = std::vector<FragmentCombination>;
        using FragmentPtr = std::shared_ptr<Fragment>;
        using FragmentList = std::vector<FragmentPtr>;

        // Constructors
        MBEFragmentGenerator() = default;
        MBEFragmentGenerator(const MBEConfig& config);

        // Main generation methods
        /**
         * @brief Generate all fragment combinations for MBE
         * @param molecule Input molecule to fragment
         * @param config MBE configuration
         * @return All fragment combinations organized by N-body order
         */
        std::unordered_map<int, FragmentCombinations> generate_all_combinations(
            const Molecule& molecule, const MBEConfig& config);

        /**
         * @brief Generate N-body fragment combinations
         * @param fragments Base fragments
         * @param n_body N-body order
         * @param config MBE configuration
         * @return All N-body combinations
         */
        FragmentCombinations generate_n_body_combinations(
            const FragmentList& fragments, int n_body, const MBEConfig& config);

        // Fragmentation methods
        /**
         * @brief Create initial fragments from molecule
         * @param molecule Input molecule
         * @param config MBE configuration
         * @return Initial fragment list
         */
        FragmentList create_initial_fragments(const Molecule& molecule, const MBEConfig& config);

        /**
         * @brief Fragment molecule using atomic fragmentation
         * @param molecule Input molecule
         * @return One fragment per atom
         */
        FragmentList fragment_by_atoms(const Molecule& molecule);

        /**
         * @brief Fragment molecule using molecular fragmentation
         * @param molecule Input molecule
         * @return One fragment per connected component
         */
        FragmentList fragment_by_molecules(const Molecule& molecule);

        /**
         * @brief Fragment molecule using functional group fragmentation
         * @param molecule Input molecule
         * @return Fragments based on functional groups
         */
        FragmentList fragment_by_functional_groups(const Molecule& molecule);

        /**
         * @brief Fragment molecule using distance-based fragmentation
         * @param molecule Input molecule
         * @param cutoff Distance cutoff for fragmentation
         * @return Distance-based fragments
         */
        FragmentList fragment_by_distance(const Molecule& molecule, double cutoff);

        /**
         * @brief Use custom fragmentation scheme
         * @param molecule Input molecule
         * @param fragment_definitions Custom fragment definitions
         * @return Custom fragments
         */
        FragmentList fragment_by_custom(const Molecule& molecule, 
            const std::vector<std::vector<std::size_t>>& fragment_definitions);

        // Filtering and truncation
        /**
         * @brief Filter combinations based on distance criteria
         * @param combinations Fragment combinations to filter
         * @param fragments Base fragments
         * @param cutoff Distance cutoff
         * @return Filtered combinations
         */
        FragmentCombinations filter_by_distance(const FragmentCombinations& combinations,
            const FragmentList& fragments, double cutoff);

        /**
         * @brief Filter combinations based on energy criteria
         * @param combinations Fragment combinations to filter
         * @param fragments Base fragments
         * @param threshold Energy threshold
         * @return Filtered combinations
         */
        FragmentCombinations filter_by_energy(const FragmentCombinations& combinations,
            const FragmentList& fragments, double threshold);

        /**
         * @brief Remove duplicate combinations
         * @param combinations Fragment combinations
         * @return Unique combinations only
         */
        FragmentCombinations remove_duplicates(const FragmentCombinations& combinations);

        // Utilities
        /**
         * @brief Calculate distance between fragment centers
         * @param frag1 First fragment
         * @param frag2 Second fragment
         * @return Distance in Angstroms
         */
        double calculate_fragment_distance(const Fragment& frag1, const Fragment& frag2);

        /**
         * @brief Calculate center of mass for a fragment combination
         * @param combination Fragment indices
         * @param fragments Base fragments
         * @return Center of mass coordinates
         */
        std::array<double, 3> calculate_combination_center(
            const FragmentCombination& combination, const FragmentList& fragments);

        /**
         * @brief Check if fragment combination is valid
         * @param combination Fragment combination to validate
         * @param fragments Base fragments
         * @param config MBE configuration
         * @return True if combination is valid
         */
        bool is_valid_combination(const FragmentCombination& combination,
            const FragmentList& fragments, const MBEConfig& config);

        /**
         * @brief Estimate number of combinations for given order
         * @param n_fragments Number of base fragments
         * @param n_body N-body order
         * @return Estimated number of combinations
         */
        std::size_t estimate_combinations(std::size_t n_fragments, int n_body);

        /**
         * @brief Get statistics about fragment generation
         * @return Map with generation statistics
         */
        std::unordered_map<std::string, std::size_t> generation_statistics() const;

        // Combination utilities
        /**
         * @brief Create fragment from combination of base fragments
         * @param combination Fragment indices
         * @param fragments Base fragments
         * @param id Fragment ID for the combination
         * @return Combined fragment
         */
        FragmentPtr create_combined_fragment(const FragmentCombination& combination,
            const FragmentList& fragments, const std::string& id);

        /**
         * @brief Generate unique ID for fragment combination
         * @param combination Fragment indices
         * @param n_body N-body order
         * @return Unique string ID
         */
        std::string generate_combination_id(const FragmentCombination& combination, int n_body);

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

    private:
        // Configuration
        MBEConfig config_;
        
        // Statistics
        mutable std::unordered_map<std::string, std::size_t> stats_;
        
        // Internal methods
        void generate_combinations_recursive(const FragmentList& fragments,
            std::size_t start_idx, int remaining_picks, 
            FragmentCombination& current_combination,
            FragmentCombinations& results);

        bool passes_distance_filter(const FragmentCombination& combination,
            const FragmentList& fragments, double cutoff);

        bool passes_energy_filter(const FragmentCombination& combination,
            const FragmentList& fragments, double threshold);

        std::vector<std::vector<std::size_t>> identify_functional_groups(const Molecule& molecule);
        std::vector<std::vector<std::size_t>> identify_connected_components(const Molecule& molecule);
        
        void update_statistics(const std::string& key, std::size_t value) const;
        
        // Functional group patterns (SMARTS-like patterns)
        struct FunctionalGroupPattern {
            std::string name;
            std::vector<int> atom_pattern;
            std::vector<std::pair<int, int>> bond_pattern;
        };
        
        std::vector<FunctionalGroupPattern> get_functional_group_patterns() const;
        bool match_functional_group(const Molecule& molecule, const FunctionalGroupPattern& pattern,
            std::size_t start_atom, std::vector<std::size_t>& matched_atoms);
    };

    /**
     * @brief Utility functions for fragment combinations
     */
    namespace mbe_utils {
        /**
         * @brief Calculate binomial coefficient (n choose k)
         * @param n Total number of items
         * @param k Number of items to choose
         * @return Binomial coefficient
         */
        std::size_t binomial_coefficient(std::size_t n, std::size_t k);

        /**
         * @brief Generate all combinations of k items from n items
         * @param n Total number of items
         * @param k Number of items to choose
         * @return All combinations as index vectors
         */
        std::vector<std::vector<std::size_t>> generate_combinations(std::size_t n, std::size_t k);
        
        /**
         * @brief Check if two fragment combinations overlap
         * @param combo1 First combination
         * @param combo2 Second combination
         * @return True if combinations share fragments
         */
        bool combinations_overlap(const std::vector<std::size_t>& combo1,
                                 const std::vector<std::size_t>& combo2);
    }

} // namespace libfrag

#endif // LIBFRAG_MBE_FRAGMENT_GENERATOR_HPP
