#pragma once
#ifndef LIBFRAG_FRAGMENT_HPP
#define LIBFRAG_FRAGMENT_HPP

#include "libfrag/molecule.hpp"
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <optional>
#include <functional>

namespace libfrag {

    // Forward declarations
    class Fragment;
    class FragmentLink;

    /**
     * @brief Represents a covalent link between two fragments
     * 
     * FragmentLink stores information about bonds that connect different fragments,
     * including the atoms involved and the bond properties.
     */
    class FragmentLink {
    public:
        using FragmentId = std::string;
        using AtomIndex = std::size_t;

        /**
         * @brief Construct a fragment link
         * @param source_fragment_id ID of the source fragment
         * @param target_fragment_id ID of the target fragment  
         * @param source_atom_index Index of atom in source fragment
         * @param target_atom_index Index of atom in target fragment
         * @param bond_type Type of bond connecting the fragments
         * @param bond_order Bond order (for fractional bonds)
         */
        FragmentLink(const FragmentId& source_fragment_id,
                    const FragmentId& target_fragment_id,
                    AtomIndex source_atom_index,
                    AtomIndex target_atom_index,
                    Bond::BondType bond_type = Bond::BondType::SINGLE,
                    double bond_order = 1.0);

        // Getters
        const FragmentId& source_fragment_id() const { return source_fragment_id_; }
        const FragmentId& target_fragment_id() const { return target_fragment_id_; }
        AtomIndex source_atom_index() const { return source_atom_index_; }
        AtomIndex target_atom_index() const { return target_atom_index_; }
        Bond::BondType bond_type() const { return bond_type_; }
        double bond_order() const { return bond_order_; }

        // Setters
        void set_bond_type(Bond::BondType bond_type) { bond_type_ = bond_type; }
        void set_bond_order(double bond_order) { bond_order_ = bond_order; }

        // Comparison operators
        bool operator==(const FragmentLink& other) const;
        bool operator!=(const FragmentLink& other) const { return !(*this == other); }

        // String representation
        std::string to_string() const;

    private:
        FragmentId source_fragment_id_;
        FragmentId target_fragment_id_;
        AtomIndex source_atom_index_;
        AtomIndex target_atom_index_;
        Bond::BondType bond_type_;
        double bond_order_;
    };

    /**
     * @brief Fragment class extending Molecule with fragment-specific functionality
     * 
     * The Fragment class implements the Composite design pattern, allowing fragments
     * to contain subfragments. It extends Molecule to add fragment identification,
     * inter-fragment connectivity, and hierarchical fragmentation capabilities.
     * 
     * Key features:
     * - Unique fragment identification system
     * - Inter-fragment covalent links management
     * - Composite pattern for hierarchical fragmentation
     * - Fragment splitting and merging operations
     * - Parent-child relationship tracking
     * - Fragment-specific properties and metadata
     */
    class Fragment : public Molecule {
    public:
        using FragmentId = std::string;
        using FragmentPtr = std::shared_ptr<Fragment>;
        using FragmentContainer = std::vector<FragmentPtr>;
        using LinkContainer = std::vector<FragmentLink>;

        // Constructors
        /**
         * @brief Default constructor creates empty fragment
         */
        Fragment() = default;

        /**
         * @brief Construct fragment with unique ID
         * @param fragment_id Unique identifier for this fragment
         */
        explicit Fragment(const FragmentId& fragment_id);

        /**
         * @brief Construct fragment from molecule with ID
         * @param molecule Molecule to copy data from
         * @param fragment_id Unique identifier for this fragment
         */
        Fragment(const Molecule& molecule, const FragmentId& fragment_id);

        /**
         * @brief Construct fragment from atoms with ID
         * @param atom_list Vector of atoms to copy
         * @param fragment_id Unique identifier for this fragment
         */
        Fragment(const std::vector<Atom>& atom_list, const FragmentId& fragment_id);

        /**
         * @brief Construct fragment with atoms, bonds, and ID
         * @param atom_list Vector of atoms to copy
         * @param bond_connectivity Vector of atom index pairs representing bonds
         * @param fragment_id Unique identifier for this fragment
         * @param bond_types Optional bond types
         */
        Fragment(const std::vector<Atom>& atom_list,
                const std::vector<std::pair<AtomIndex, AtomIndex>>& bond_connectivity,
                const FragmentId& fragment_id,
                const std::vector<Bond::BondType>& bond_types = {});

        // Copy and move constructors
        Fragment(const Fragment& other);
        Fragment& operator=(const Fragment& other);
        Fragment(Fragment&& other) = default;
        Fragment& operator=(Fragment&& other) = default;

        // Destructor
        virtual ~Fragment() = default;

        // Fragment identification
        /**
         * @brief Get fragment ID
         * @return Unique fragment identifier
         */
        const FragmentId& fragment_id() const { return fragment_id_; }

        /**
         * @brief Set fragment ID
         * @param new_id New unique identifier
         */
        void set_fragment_id(const FragmentId& new_id) { fragment_id_ = new_id; }

        /**
         * @brief Get fragment generation level (0 for root, 1 for first split, etc.)
         * @return Generation level in fragment hierarchy
         */
        std::size_t generation_level() const { return generation_level_; }

        // Parent-child relationship management
        /**
         * @brief Get parent fragment
         * @return Weak pointer to parent fragment (empty if root)
         */
        std::weak_ptr<Fragment> parent_fragment() const { return parent_fragment_; }

        /**
         * @brief Set parent fragment
         * @param parent Shared pointer to parent fragment
         */
        void set_parent_fragment(FragmentPtr parent);

        /**
         * @brief Check if this is a root fragment (no parent)
         * @return True if fragment has no parent
         */
        bool is_root_fragment() const { return parent_fragment_.expired(); }

        /**
         * @brief Get root fragment in hierarchy
         * @return Shared pointer to root fragment
         */
        FragmentPtr get_root_fragment();

        // Subfragment management (Composite pattern)
        /**
         * @brief Get number of direct subfragments
         * @return Number of child fragments
         */
        std::size_t subfragment_count() const { return subfragments_.size(); }

        /**
         * @brief Get all direct subfragments
         * @return Vector of subfragment pointers
         */
        const FragmentContainer& subfragments() const { return subfragments_; }

        /**
         * @brief Get subfragment by index
         * @param index Index of subfragment
         * @return Shared pointer to subfragment
         * @throws std::out_of_range if index is invalid
         */
        FragmentPtr subfragment(std::size_t index) const;

        /**
         * @brief Add subfragment to this fragment
         * @param subfragment Shared pointer to subfragment to add
         */
        void add_subfragment(FragmentPtr subfragment);

        /**
         * @brief Remove subfragment by index
         * @param index Index of subfragment to remove
         * @throws std::out_of_range if index is invalid
         */
        void remove_subfragment(std::size_t index);

        /**
         * @brief Remove all subfragments
         */
        void clear_subfragments();

        /**
         * @brief Check if fragment has subfragments
         * @return True if fragment contains subfragments
         */
        bool has_subfragments() const { return !subfragments_.empty(); }

        /**
         * @brief Get all fragments in subtree (recursive)
         * @return Vector of all fragments in hierarchy below this one
         */
        FragmentContainer get_all_subfragments() const;

        /**
         * @brief Get total number of atoms in entire fragment tree
         * @return Total atom count including all subfragments
         */
        std::size_t total_atom_count() const;

        // Inter-fragment connectivity
        /**
         * @brief Get all fragment links
         * @return Vector of fragment links
         */
        const LinkContainer& fragment_links() const { return fragment_links_; }

        /**
         * @brief Add link to another fragment
         * @param link FragmentLink to add
         */
        void add_fragment_link(const FragmentLink& link);

        /**
         * @brief Remove fragment link by index
         * @param index Index of link to remove
         * @throws std::out_of_range if index is invalid
         */
        void remove_fragment_link(std::size_t index);

        /**
         * @brief Find links to specific fragment
         * @param target_fragment_id ID of target fragment
         * @return Vector of indices of links to target fragment
         */
        std::vector<std::size_t> find_links_to_fragment(const FragmentId& target_fragment_id) const;

        /**
         * @brief Get all fragment IDs this fragment is linked to
         * @return Set of fragment IDs
         */
        std::unordered_set<FragmentId> linked_fragment_ids() const;

        /**
         * @brief Check if fragment is linked to another fragment
         * @param target_fragment_id ID of target fragment
         * @return True if fragments are linked
         */
        bool is_linked_to(const FragmentId& target_fragment_id) const;

        // Fragment splitting operations
        /**
         * @brief Split fragment into subfragments based on connectivity
         * @param max_fragment_size Maximum number of atoms per subfragment (0 = no limit)
         * @return Vector of newly created subfragments
         */
        FragmentContainer fragment_by_connectivity(std::size_t max_fragment_size = 0);

        /**
         * @brief Split fragment by breaking specific bonds
         * @param bonds_to_break Indices of bonds to break for fragmentation
         * @return Vector of newly created subfragments
         */
        FragmentContainer fragment_by_bond_breaking(const std::vector<BondIndex>& bonds_to_break);

        /**
         * @brief Split fragment using custom fragmentation function
         * @param fragmenter Function that takes atoms and bonds and returns fragment assignments
         * @return Vector of newly created subfragments
         */
        FragmentContainer fragment_by_custom_function(
            const std::function<std::vector<std::vector<AtomIndex>>(const Fragment&)>& fragmenter);

        /**
         * @brief Split fragment into roughly equal-sized pieces
         * @param target_fragment_count Desired number of subfragments
         * @return Vector of newly created subfragments
         */
        FragmentContainer fragment_by_size(std::size_t target_fragment_count);

        // Fragment merging operations
        /**
         * @brief Merge all subfragments back into this fragment
         * @param preserve_links Whether to preserve inter-fragment links as internal bonds
         */
        void merge_subfragments(bool preserve_links = true);

        /**
         * @brief Merge specific subfragments
         * @param subfragment_indices Indices of subfragments to merge
         * @param preserve_links Whether to preserve inter-fragment links as internal bonds
         */
        void merge_selected_subfragments(const std::vector<std::size_t>& subfragment_indices,
                                       bool preserve_links = true);

        // Fragment analysis
        /**
         * @brief Calculate fragment-specific properties
         * @return Map of fragment property names to values
         */
        std::unordered_map<std::string, double> calculate_fragment_properties() const;

        /**
         * @brief Get fragment composition (element counts)
         * @return Map of element symbols to counts
         */
        std::unordered_map<std::string, std::size_t> fragment_composition() const;

        /**
         * @brief Generate fragment hash for quick comparison
         * @return Hash string based on composition and connectivity
         */
        std::string fragment_hash() const;

        // String representations
        /**
         * @brief Generate detailed fragment information
         * @return Multi-line string with fragment hierarchy and properties
         */
        std::string to_fragment_string() const;

        /**
         * @brief Generate fragment tree representation
         * @return Multi-line string showing fragment hierarchy
         */
        std::string to_tree_string(std::size_t indent_level = 0) const;

        // Static factory methods
        /**
         * @brief Create fragment from molecule with automatic ID generation
         * @param molecule Source molecule
         * @param id_prefix Prefix for generated ID
         * @return New fragment with generated ID
         */
        static FragmentPtr create_from_molecule(const Molecule& molecule, 
                                              const std::string& id_prefix = "frag");

        /**
         * @brief Create fragment hierarchy from molecule using automatic fragmentation
         * @param molecule Source molecule
         * @param max_fragment_size Maximum atoms per fragment
         * @param id_prefix Prefix for generated IDs
         * @return Root fragment with subfragments
         */
        static FragmentPtr create_fragment_hierarchy(const Molecule& molecule,
                                                   std::size_t max_fragment_size = 10,
                                                   const std::string& id_prefix = "frag");

        // Operators
        bool operator==(const Fragment& other) const;
        bool operator!=(const Fragment& other) const { return !(*this == other); }

    private:
        // Fragment-specific data
        FragmentId fragment_id_;                        ///< Unique fragment identifier
        std::size_t generation_level_ = 0;              ///< Level in fragment hierarchy
        std::weak_ptr<Fragment> parent_fragment_;       ///< Parent fragment (empty for root)
        FragmentContainer subfragments_;                ///< Direct child fragments
        LinkContainer fragment_links_;                  ///< Links to other fragments
        
        // Fragment metadata
        std::unordered_map<std::string, std::string> fragment_metadata_;  ///< String metadata
        
        // Static ID generation
        static std::size_t next_fragment_id_;
        
        // Helper methods
        void validate_subfragment_index(std::size_t index) const;
        FragmentId generate_child_id(const std::string& suffix = "") const;
        void update_generation_levels();
        void collect_all_subfragments(FragmentContainer& all_fragments) const;
        
        // Fragmentation helpers
        std::vector<std::vector<AtomIndex>> find_connected_components() const;
        FragmentPtr create_subfragment_from_atoms(const std::vector<AtomIndex>& atom_indices,
                                                 const std::string& suffix) const;
        void transfer_relevant_bonds(const std::vector<AtomIndex>& atom_indices,
                                   Fragment& target_fragment,
                                   std::unordered_map<AtomIndex, AtomIndex>& atom_index_map) const;
    };

    // Inline method implementations

    inline FragmentLink::FragmentLink(const FragmentId& source_fragment_id,
                                     const FragmentId& target_fragment_id,
                                     AtomIndex source_atom_index,
                                     AtomIndex target_atom_index,
                                     Bond::BondType bond_type,
                                     double bond_order)
        : source_fragment_id_(source_fragment_id)
        , target_fragment_id_(target_fragment_id)
        , source_atom_index_(source_atom_index)
        , target_atom_index_(target_atom_index)
        , bond_type_(bond_type)
        , bond_order_(bond_order) {}

    inline Fragment::FragmentPtr Fragment::subfragment(std::size_t index) const {
        validate_subfragment_index(index);
        return subfragments_[index];
    }

    inline void Fragment::validate_subfragment_index(std::size_t index) const {
        if (index >= subfragments_.size()) {
            throw std::out_of_range("Subfragment index " + std::to_string(index) + 
                                  " out of range (fragment has " + std::to_string(subfragments_.size()) + " subfragments)");
        }
    }

    inline std::unordered_set<Fragment::FragmentId> Fragment::linked_fragment_ids() const {
        std::unordered_set<FragmentId> linked_ids;
        for (const auto& link : fragment_links_) {
            linked_ids.insert(link.target_fragment_id());
        }
        return linked_ids;
    }

    inline bool Fragment::is_linked_to(const FragmentId& target_fragment_id) const {
        return !find_links_to_fragment(target_fragment_id).empty();
    }

} // namespace libfrag

#endif // LIBFRAG_FRAGMENT_HPP
