#pragma once
#ifndef LIBFRAG_GLOBAL_SCF_FRAGMENT_GENERATOR_HPP
#define LIBFRAG_GLOBAL_SCF_FRAGMENT_GENERATOR_HPP

#include <vector>
#include <memory>
#include <unordered_map>
#include <string>
#include "molecule.hpp"
#include "fragment.hpp"
#include "global_scf_config.hpp"

namespace libfrag {

    /**
     * @brief Fragment definition for Global SCF calculations
     * 
     * Extended fragment information including buffer regions and
     * embedding environment for Global SCF calculations.
     */
    struct GlobalSCFFragment {
        // Core fragment data
        std::shared_ptr<Fragment> core_fragment;   ///< Core fragment object
        std::size_t fragment_id;                   ///< Unique fragment identifier
        std::string fragment_name;                 ///< Human-readable name
        
        // Atom assignments
        std::vector<std::size_t> core_atoms;       ///< Core fragment atoms
        std::vector<std::size_t> buffer_atoms;     ///< Buffer region atoms
        std::vector<std::size_t> environment_atoms; ///< Environment atoms for embedding
        
        // Geometric properties
        std::array<double, 3> center_of_mass;      ///< Fragment center of mass
        double radius;                             ///< Fragment radius (approximate)
        std::array<double, 6> bounding_box;        ///< Fragment bounding box
        
        // Fragment relationships
        std::vector<std::size_t> neighboring_fragments; ///< IDs of neighboring fragments
        std::unordered_map<std::size_t, double> fragment_distances; ///< Distances to other fragments
        
        // Buffer and embedding information
        double buffer_distance;                    ///< Buffer region distance
        std::vector<std::size_t> capping_atoms;    ///< Atoms for capping dangling bonds
        
        // Charge and multiplicity
        int total_charge = 0;                      ///< Total fragment charge
        int multiplicity = 1;                      ///< Fragment spin multiplicity
        
        // Constructor
        GlobalSCFFragment() = default;
        GlobalSCFFragment(std::size_t id, std::shared_ptr<Fragment> fragment);
        
        // Utilities
        std::string to_string() const;
        bool is_valid() const;
        std::size_t total_atoms() const;
        bool contains_atom(std::size_t atom_index) const;
        bool is_neighbor(std::size_t other_fragment_id) const;
        double distance_to_fragment(std::size_t other_fragment_id) const;
    };

    /**
     * @brief Generates fragment definitions for Global SCF calculations
     * 
     * This class is responsible for creating fragment definitions suitable for
     * Global SCF calculations, including core fragments, buffer regions, and
     * embedding environments. It handles various fragmentation schemes and
     * ensures proper fragment boundaries and relationships.
     */
    class GlobalSCFFragmentGenerator {
    public:
        // Type aliases
        using FragmentList = std::vector<GlobalSCFFragment>;
        using NeighborMap = std::unordered_map<std::size_t, std::vector<std::size_t>>;
        using DistanceMatrix = std::unordered_map<std::pair<std::size_t, std::size_t>, double>;

        // Constructors
        GlobalSCFFragmentGenerator() = default;
        explicit GlobalSCFFragmentGenerator(const GlobalSCFConfig& config);

        // Configuration
        /**
         * @brief Set configuration
         * @param config Global SCF configuration
         */
        void set_config(const GlobalSCFConfig& config) { config_ = config; }

        /**
         * @brief Get current configuration
         * @return Current configuration
         */
        const GlobalSCFConfig& config() const { return config_; }

        // Main fragmentation methods
        /**
         * @brief Generate fragment definitions for molecule
         * @param molecule Input molecule
         * @return List of fragment definitions
         */
        FragmentList generate_fragments(const Molecule& molecule);

        /**
         * @brief Generate fragments with specific configuration
         * @param molecule Input molecule
         * @param config Global SCF configuration
         * @return List of fragment definitions
         */
        FragmentList generate_fragments(const Molecule& molecule, const GlobalSCFConfig& config);

        // Specific fragmentation schemes
        /**
         * @brief Fragment by individual atoms
         * @param molecule Input molecule
         * @return List of atomic fragments
         */
        FragmentList fragment_by_atoms(const Molecule& molecule);

        /**
         * @brief Fragment by molecular units
         * @param molecule Input molecule
         * @return List of molecular fragments
         */
        FragmentList fragment_by_molecules(const Molecule& molecule);

        /**
         * @brief Fragment by functional groups
         * @param molecule Input molecule
         * @return List of functional group fragments
         */
        FragmentList fragment_by_functional_groups(const Molecule& molecule);

        /**
         * @brief Fragment using distance-based clustering
         * @param molecule Input molecule
         * @param cutoff Distance cutoff for clustering
         * @return List of distance-based fragments
         */
        FragmentList fragment_by_distance(const Molecule& molecule, double cutoff);

        /**
         * @brief Fragment using custom definitions
         * @param molecule Input molecule
         * @param fragment_definitions User-defined fragment atom lists
         * @return List of custom fragments
         */
        FragmentList fragment_by_custom(const Molecule& molecule, 
                                       const std::vector<std::vector<std::size_t>>& fragment_definitions);

        /**
         * @brief Fragment using FMO-style approach
         * @param molecule Input molecule
         * @return List of FMO-style fragments
         */
        FragmentList fragment_fmo_style(const Molecule& molecule);

        // Buffer and environment assignment
        /**
         * @brief Assign buffer regions to fragments
         * @param fragments Fragment list to modify
         * @param molecule Source molecule
         * @param buffer_distance Buffer distance in Angstroms
         */
        void assign_buffer_regions(FragmentList& fragments, const Molecule& molecule, 
                                  double buffer_distance);

        /**
         * @brief Assign environment atoms for embedding
         * @param fragments Fragment list to modify
         * @param molecule Source molecule
         */
        void assign_environment_atoms(FragmentList& fragments, const Molecule& molecule);

        /**
         * @brief Add capping atoms for dangling bonds
         * @param fragments Fragment list to modify
         * @param molecule Source molecule
         */
        void add_capping_atoms(FragmentList& fragments, const Molecule& molecule);

        // Fragment relationship analysis
        /**
         * @brief Determine fragment neighbors
         * @param fragments Fragment list
         * @param molecule Source molecule
         * @return Neighbor mapping
         */
        NeighborMap determine_neighbors(const FragmentList& fragments, const Molecule& molecule);

        /**
         * @brief Calculate fragment-fragment distances
         * @param fragments Fragment list
         * @return Distance matrix
         */
        std::unordered_map<std::string, double> calculate_distances(const FragmentList& fragments);

        /**
         * @brief Analyze fragment overlap
         * @param fragments Fragment list
         * @return Overlap analysis results
         */
        std::unordered_map<std::string, double> analyze_overlap(const FragmentList& fragments);

        // Validation and optimization
        /**
         * @brief Validate fragment definitions
         * @param fragments Fragment list
         * @param molecule Source molecule
         * @return List of validation issues
         */
        std::vector<std::string> validate_fragments(const FragmentList& fragments, 
                                                    const Molecule& molecule);

        /**
         * @brief Optimize fragment boundaries
         * @param fragments Fragment list to optimize
         * @param molecule Source molecule
         */
        void optimize_boundaries(FragmentList& fragments, const Molecule& molecule);

        /**
         * @brief Balance fragment sizes
         * @param fragments Fragment list to balance
         * @param molecule Source molecule
         * @param target_size Target fragment size
         */
        void balance_fragment_sizes(FragmentList& fragments, const Molecule& molecule, 
                                   std::size_t target_size);

        // Statistics and analysis
        /**
         * @brief Get fragmentation statistics
         * @param fragments Fragment list
         * @return Statistics map
         */
        std::unordered_map<std::string, double> get_fragmentation_statistics(
            const FragmentList& fragments) const;

        /**
         * @brief Estimate computational cost
         * @param fragments Fragment list
         * @param config Global SCF configuration
         * @return Cost estimate map
         */
        std::unordered_map<std::string, double> estimate_computational_cost(
            const FragmentList& fragments, const GlobalSCFConfig& config) const;

        /**
         * @brief Generate fragmentation report
         * @param fragments Fragment list
         * @param molecule Source molecule
         * @return Formatted report string
         */
        std::string generate_report(const FragmentList& fragments, const Molecule& molecule) const;

        // Utility methods
        /**
         * @brief Create fragment from atom indices
         * @param molecule Source molecule
         * @param atom_indices Atom indices for fragment
         * @param fragment_id Fragment identifier
         * @return GlobalSCFFragment object
         */
        GlobalSCFFragment create_fragment_from_atoms(const Molecule& molecule,
                                                    const std::vector<std::size_t>& atom_indices,
                                                    std::size_t fragment_id);

        /**
         * @brief Find atoms within distance of fragment
         * @param molecule Source molecule
         * @param fragment Fragment to search around
         * @param distance Search distance
         * @return List of atom indices
         */
        std::vector<std::size_t> find_atoms_within_distance(const Molecule& molecule,
                                                           const GlobalSCFFragment& fragment,
                                                           double distance);

        /**
         * @brief Get fragment containing atom
         * @param fragments Fragment list
         * @param atom_index Atom index to find
         * @return Fragment ID (or SIZE_MAX if not found)
         */
        std::size_t get_fragment_containing_atom(const FragmentList& fragments, 
                                                std::size_t atom_index);

    private:
        // Configuration
        GlobalSCFConfig config_;
        
        // Internal state
        mutable std::unordered_map<std::string, double> statistics_;
        
        // Helper methods for fragmentation schemes
        std::vector<std::vector<std::size_t>> identify_molecular_units(const Molecule& molecule);
        std::vector<std::vector<std::size_t>> identify_functional_groups(const Molecule& molecule);
        std::vector<std::vector<std::size_t>> cluster_by_distance(const Molecule& molecule, double cutoff);
        std::vector<std::vector<std::size_t>> create_fmo_fragments(const Molecule& molecule);
        
        // Helper methods for buffer and environment
        std::vector<std::size_t> find_buffer_atoms(const Molecule& molecule,
                                                  const std::vector<std::size_t>& core_atoms,
                                                  double buffer_distance);
        std::vector<std::size_t> find_environment_atoms(const Molecule& molecule,
                                                       const std::vector<std::size_t>& core_atoms,
                                                       const std::vector<std::size_t>& buffer_atoms);
        
        // Helper methods for capping
        std::vector<std::size_t> identify_dangling_bonds(const Molecule& molecule,
                                                        const std::vector<std::size_t>& fragment_atoms);
        void add_hydrogen_caps(const Molecule& molecule, GlobalSCFFragment& fragment);
        
        // Helper methods for validation
        bool check_fragment_completeness(const FragmentList& fragments, const Molecule& molecule);
        bool check_fragment_overlap(const FragmentList& fragments);
        bool check_fragment_connectivity(const GlobalSCFFragment& fragment, const Molecule& molecule);
        
        // Utility functions
        double calculate_fragment_distance(const GlobalSCFFragment& frag1, const GlobalSCFFragment& frag2);
        std::array<double, 3> calculate_fragment_center(const Molecule& molecule,
                                                       const std::vector<std::size_t>& atom_indices);
        double calculate_fragment_radius(const Molecule& molecule,
                                        const std::vector<std::size_t>& atom_indices,
                                        const std::array<double, 3>& center);
        std::array<double, 6> calculate_bounding_box(const Molecule& molecule,
                                                    const std::vector<std::size_t>& atom_indices);
        
        // Statistical helpers
        void update_statistics(const std::string& key, double value) const;
        void increment_statistics(const std::string& key, double increment = 1.0) const;
    };

}

#endif // LIBFRAG_GLOBAL_SCF_FRAGMENT_GENERATOR_HPP
