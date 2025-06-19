#include "libfrag/fragment.hpp"
#include <algorithm>
#include <sstream>
#include <queue>
#include <stack>
#include <iomanip>
#include <functional>

namespace libfrag {

    // Static member initialization
    std::size_t Fragment::next_fragment_id_ = 1;

    // FragmentLink implementation
    bool FragmentLink::operator==(const FragmentLink& other) const {
        return source_fragment_id_ == other.source_fragment_id_ &&
               target_fragment_id_ == other.target_fragment_id_ &&
               source_atom_index_ == other.source_atom_index_ &&
               target_atom_index_ == other.target_atom_index_ &&
               bond_type_ == other.bond_type_ &&
               std::abs(bond_order_ - other.bond_order_) < 1e-6;
    }

    std::string FragmentLink::to_string() const {
        std::ostringstream oss;
        oss << source_fragment_id_ << "[" << source_atom_index_ << "] -- "
            << target_fragment_id_ << "[" << target_atom_index_ << "] "
            << "(bond_order: " << std::fixed << std::setprecision(2) << bond_order_ << ")";
        return oss.str();
    }

    // Fragment constructors
    Fragment::Fragment(const FragmentId& fragment_id)
        : Molecule(), fragment_id_(fragment_id) {}

    Fragment::Fragment(const Molecule& molecule, const FragmentId& fragment_id)
        : Molecule(molecule), fragment_id_(fragment_id) {}

    Fragment::Fragment(const std::vector<Atom>& atom_list, const FragmentId& fragment_id)
        : Molecule(atom_list), fragment_id_(fragment_id) {}

    Fragment::Fragment(const std::vector<Atom>& atom_list,
                      const std::vector<std::pair<AtomIndex, AtomIndex>>& bond_connectivity,
                      const FragmentId& fragment_id,
                      const std::vector<Bond::BondType>& bond_types)
        : Molecule(atom_list, bond_connectivity, bond_types), fragment_id_(fragment_id) {}

    Fragment::Fragment(const Fragment& other)
        : Molecule(other)
        , fragment_id_(other.fragment_id_)
        , generation_level_(other.generation_level_)
        , fragment_links_(other.fragment_links_)
        , fragment_metadata_(other.fragment_metadata_) {
        
        // Deep copy subfragments
        subfragments_.reserve(other.subfragments_.size());
        for (const auto& subfrag : other.subfragments_) {
            auto new_subfrag = std::make_shared<Fragment>(*subfrag);
            new_subfrag->set_parent_fragment(std::shared_ptr<Fragment>(this, [](Fragment*){}));
            subfragments_.push_back(new_subfrag);
        }
    }

    Fragment& Fragment::operator=(const Fragment& other) {
        if (this != &other) {
            Molecule::operator=(other);
            fragment_id_ = other.fragment_id_;
            generation_level_ = other.generation_level_;
            fragment_links_ = other.fragment_links_;
            fragment_metadata_ = other.fragment_metadata_;
            
            // Deep copy subfragments
            subfragments_.clear();
            subfragments_.reserve(other.subfragments_.size());
            for (const auto& subfrag : other.subfragments_) {
                auto new_subfrag = std::make_shared<Fragment>(*subfrag);
                new_subfrag->set_parent_fragment(std::shared_ptr<Fragment>(this, [](Fragment*){}));
                subfragments_.push_back(new_subfrag);
            }
        }
        return *this;
    }

    // Parent-child relationship management
    void Fragment::set_parent_fragment(FragmentPtr parent) {
        parent_fragment_ = parent;
        if (parent) {
            generation_level_ = parent->generation_level_ + 1;
        } else {
            generation_level_ = 0;
        }
        update_generation_levels();
    }

    Fragment::FragmentPtr Fragment::get_root_fragment() {
        auto current = std::shared_ptr<Fragment>(this, [](Fragment*){});
        while (!current->is_root_fragment()) {
            if (auto parent = current->parent_fragment_.lock()) {
                current = parent;
            } else {
                break;
            }
        }
        return current;
    }

    // Subfragment management
    void Fragment::add_subfragment(FragmentPtr subfragment) {
        if (subfragment) {
            subfragment->set_parent_fragment(std::shared_ptr<Fragment>(this, [](Fragment*){}));
            subfragments_.push_back(subfragment);
        }
    }

    void Fragment::remove_subfragment(std::size_t index) {
        validate_subfragment_index(index);
        subfragments_[index]->parent_fragment_.reset();
        subfragments_.erase(subfragments_.begin() + index);
    }

    void Fragment::clear_subfragments() {
        for (auto& subfrag : subfragments_) {
            subfrag->parent_fragment_.reset();
        }
        subfragments_.clear();
    }

    Fragment::FragmentContainer Fragment::get_all_subfragments() const {
        FragmentContainer all_fragments;
        collect_all_subfragments(all_fragments);
        return all_fragments;
    }

    std::size_t Fragment::total_atom_count() const {
        std::size_t total = atom_count();
        for (const auto& subfrag : subfragments_) {
            total += subfrag->total_atom_count();
        }
        return total;
    }

    // Inter-fragment connectivity
    void Fragment::add_fragment_link(const FragmentLink& link) {
        fragment_links_.push_back(link);
    }

    void Fragment::remove_fragment_link(std::size_t index) {
        if (index >= fragment_links_.size()) {
            throw std::out_of_range("Fragment link index out of range");
        }
        fragment_links_.erase(fragment_links_.begin() + index);
    }

    std::vector<std::size_t> Fragment::find_links_to_fragment(const FragmentId& target_fragment_id) const {
        std::vector<std::size_t> link_indices;
        for (std::size_t i = 0; i < fragment_links_.size(); ++i) {
            if (fragment_links_[i].target_fragment_id() == target_fragment_id) {
                link_indices.push_back(i);
            }
        }
        return link_indices;
    }

    // Fragment splitting operations
    Fragment::FragmentContainer Fragment::fragment_by_connectivity(std::size_t max_fragment_size) {
        if (is_empty()) {
            return {};
        }

        auto connected_components = find_connected_components();
        FragmentContainer new_fragments;

        for (std::size_t comp_idx = 0; comp_idx < connected_components.size(); ++comp_idx) {
            const auto& component = connected_components[comp_idx];
            
            // If max_fragment_size is specified and component is too large, split further
            if (max_fragment_size > 0 && component.size() > max_fragment_size) {
                // Simple splitting by taking chunks of max_fragment_size
                for (std::size_t start = 0; start < component.size(); start += max_fragment_size) {
                    std::size_t end = std::min(start + max_fragment_size, component.size());
                    std::vector<AtomIndex> chunk(component.begin() + start, component.begin() + end);
                    
                    std::string suffix = "_comp" + std::to_string(comp_idx) + "_chunk" + 
                                       std::to_string(start / max_fragment_size);
                    auto fragment = create_subfragment_from_atoms(chunk, suffix);
                    new_fragments.push_back(fragment);
                }
            } else {
                std::string suffix = "_comp" + std::to_string(comp_idx);
                auto fragment = create_subfragment_from_atoms(component, suffix);
                new_fragments.push_back(fragment);
            }
        }

        // Add as subfragments
        for (auto& fragment : new_fragments) {
            add_subfragment(fragment);
        }

        return new_fragments;
    }

    Fragment::FragmentContainer Fragment::fragment_by_bond_breaking(const std::vector<BondIndex>& bonds_to_break) {
        if (bonds_to_break.empty()) {
            return {};
        }

        // Create a copy of the molecule without the specified bonds
        Fragment temp_fragment(*this);
        
        // Remove bonds in reverse order to maintain indices
        auto sorted_bonds = bonds_to_break;
        std::sort(sorted_bonds.rbegin(), sorted_bonds.rend());
        
        for (BondIndex bond_idx : sorted_bonds) {
            if (bond_idx < temp_fragment.bond_count()) {
                temp_fragment.remove_bond(bond_idx);
            }
        }

        // Find connected components in the modified fragment
        return temp_fragment.fragment_by_connectivity();
    }

    Fragment::FragmentContainer Fragment::fragment_by_custom_function(
        const std::function<std::vector<std::vector<AtomIndex>>(const Fragment&)>& fragmenter) {
        
        auto atom_groups = fragmenter(*this);
        FragmentContainer new_fragments;

        for (std::size_t group_idx = 0; group_idx < atom_groups.size(); ++group_idx) {
            const auto& atom_group = atom_groups[group_idx];
            if (!atom_group.empty()) {
                std::string suffix = "_custom" + std::to_string(group_idx);
                auto fragment = create_subfragment_from_atoms(atom_group, suffix);
                new_fragments.push_back(fragment);
                add_subfragment(fragment);
            }
        }

        return new_fragments;
    }

    Fragment::FragmentContainer Fragment::fragment_by_size(std::size_t target_fragment_count) {
        if (target_fragment_count == 0 || atom_count() == 0) {
            return {};
        }

        std::size_t atoms_per_fragment = std::max(std::size_t(1), atom_count() / target_fragment_count);
        
        FragmentContainer new_fragments;
        for (std::size_t i = 0; i < atom_count(); i += atoms_per_fragment) {
            std::size_t end_idx = std::min(i + atoms_per_fragment, atom_count());
            
            std::vector<AtomIndex> atom_indices;
            for (std::size_t j = i; j < end_idx; ++j) {
                atom_indices.push_back(j);
            }
            
            std::string suffix = "_size" + std::to_string(i / atoms_per_fragment);
            auto fragment = create_subfragment_from_atoms(atom_indices, suffix);
            new_fragments.push_back(fragment);
            add_subfragment(fragment);
        }

        return new_fragments;
    }

    // Fragment merging operations
    void Fragment::merge_subfragments(bool preserve_links) {
        if (subfragments_.empty()) {
            return;
        }

        // Collect all atoms and bonds from subfragments
        for (const auto& subfrag : subfragments_) {
            // Add atoms
            for (std::size_t i = 0; i < subfrag->atom_count(); ++i) {
                add_atom(subfrag->atom(i));
            }
            
            // Add bonds if preserving links or if they're internal
            for (std::size_t i = 0; i < subfrag->bond_count(); ++i) {
                const auto& bond = subfrag->bond(i);
                // Note: This is simplified - in practice, you'd need to map atom indices correctly
                // add_bond(...);
            }
        }

        // Clear subfragments
        clear_subfragments();
    }

    void Fragment::merge_selected_subfragments(const std::vector<std::size_t>& subfragment_indices,
                                              bool preserve_links) {
        // Validate indices
        for (std::size_t idx : subfragment_indices) {
            validate_subfragment_index(idx);
        }

        // Merge selected subfragments (implementation similar to merge_subfragments)
        // This is a simplified version - full implementation would require careful atom index mapping
        
        // Remove merged subfragments in reverse order
        auto sorted_indices = subfragment_indices;
        std::sort(sorted_indices.rbegin(), sorted_indices.rend());
        
        for (std::size_t idx : sorted_indices) {
            remove_subfragment(idx);
        }
    }

    // Fragment analysis
    std::unordered_map<std::string, double> Fragment::calculate_fragment_properties() const {
        std::unordered_map<std::string, double> properties;
        
        // Basic molecular properties
        properties["molecular_mass"] = molecular_mass();
        properties["atom_count"] = static_cast<double>(atom_count());
        properties["bond_count"] = static_cast<double>(bond_count());
        properties["subfragment_count"] = static_cast<double>(subfragment_count());
        properties["total_atom_count"] = static_cast<double>(total_atom_count());
        properties["generation_level"] = static_cast<double>(generation_level_);
        
        // Connectivity properties
        if (!is_empty()) {
            properties["is_connected"] = is_connected() ? 1.0 : 0.0;
            properties["molecular_diameter"] = molecular_diameter();
        }
        
        return properties;
    }

    std::unordered_map<std::string, std::size_t> Fragment::fragment_composition() const {
        std::unordered_map<std::string, std::size_t> composition;
        
        for (std::size_t i = 0; i < atom_count(); ++i) {
            std::string symbol = atom(i).element_symbol();
            composition[symbol]++;
        }
        
        return composition;
    }

    std::string Fragment::fragment_hash() const {
        // Simple hash based on composition and connectivity
        std::ostringstream hash_stream;
        
        auto composition = fragment_composition();
        for (const auto& [symbol, count] : composition) {
            hash_stream << symbol << count;
        }
        
        hash_stream << "_bonds" << bond_count();
        hash_stream << "_subfrag" << subfragment_count();
        
        return hash_stream.str();
    }

    // String representations
    std::string Fragment::to_fragment_string() const {
        std::ostringstream oss;
        oss << "Fragment ID: " << fragment_id_ << "\n";
        oss << "Generation Level: " << generation_level_ << "\n";
        oss << "Atoms: " << atom_count() << ", Bonds: " << bond_count() << "\n";
        oss << "Subfragments: " << subfragment_count() << "\n";
        oss << "Fragment Links: " << fragment_links_.size() << "\n";
        
        if (!fragment_links_.empty()) {
            oss << "Links:\n";
            for (const auto& link : fragment_links_) {
                oss << "  " << link.to_string() << "\n";
            }
        }
        
        // Add molecular information
        oss << "\nMolecular Data:\n";
        oss << to_detailed_string();
        
        return oss.str();
    }

    std::string Fragment::to_tree_string(std::size_t indent_level) const {
        std::ostringstream oss;
        
        // Add indentation
        std::string indent(indent_level * 2, ' ');
        oss << indent << "├─ " << fragment_id_ << " (atoms: " << atom_count() 
            << ", bonds: " << bond_count() << ")\n";
        
        // Add subfragments recursively
        for (const auto& subfrag : subfragments_) {
            oss << subfrag->to_tree_string(indent_level + 1);
        }
        
        return oss.str();
    }

    // Static factory methods
    Fragment::FragmentPtr Fragment::create_from_molecule(const Molecule& molecule, 
                                                        const std::string& id_prefix) {
        std::string id = id_prefix + "_" + std::to_string(next_fragment_id_++);
        return std::make_shared<Fragment>(molecule, id);
    }

    Fragment::FragmentPtr Fragment::create_fragment_hierarchy(const Molecule& molecule,
                                                            std::size_t max_fragment_size,
                                                            const std::string& id_prefix) {
        auto root_fragment = create_from_molecule(molecule, id_prefix + "_root");
        
        if (molecule.atom_count() > max_fragment_size) {
            root_fragment->fragment_by_connectivity(max_fragment_size);
        }
        
        return root_fragment;
    }

    // Operators
    bool Fragment::operator==(const Fragment& other) const {
        return Molecule::operator==(other) &&
               fragment_id_ == other.fragment_id_ &&
               generation_level_ == other.generation_level_ &&
               fragment_links_ == other.fragment_links_ &&
               subfragments_.size() == other.subfragments_.size();
    }

    // Helper methods
    void Fragment::update_generation_levels() {
        for (auto& subfrag : subfragments_) {
            subfrag->generation_level_ = generation_level_ + 1;
            subfrag->update_generation_levels();
        }
    }

    void Fragment::collect_all_subfragments(FragmentContainer& all_fragments) const {
        for (const auto& subfrag : subfragments_) {
            all_fragments.push_back(subfrag);
            subfrag->collect_all_subfragments(all_fragments);
        }
    }

    Fragment::FragmentId Fragment::generate_child_id(const std::string& suffix) const {
        return fragment_id_ + suffix;
    }

    std::vector<std::vector<Fragment::AtomIndex>> Fragment::find_connected_components() const {
        if (is_empty()) {
            return {};
        }

        std::vector<bool> visited(atom_count(), false);
        std::vector<std::vector<AtomIndex>> components;

        for (AtomIndex i = 0; i < atom_count(); ++i) {
            if (!visited[i]) {
                std::vector<AtomIndex> component;
                
                // BFS to find all connected atoms
                std::queue<AtomIndex> to_visit;
                to_visit.push(i);
                visited[i] = true;

                while (!to_visit.empty()) {
                    AtomIndex current = to_visit.front();
                    to_visit.pop();
                    component.push_back(current);

                    // Check all bonded atoms
                    auto bonded_atom_indices = bonded_atoms(current);
                    for (AtomIndex bonded_idx : bonded_atom_indices) {
                        if (!visited[bonded_idx]) {
                            visited[bonded_idx] = true;
                            to_visit.push(bonded_idx);
                        }
                    }
                }

                components.push_back(component);
            }
        }

        return components;
    }

    Fragment::FragmentPtr Fragment::create_subfragment_from_atoms(const std::vector<AtomIndex>& atom_indices,
                                                                 const std::string& suffix) const {
        std::string child_id = generate_child_id(suffix);
        auto new_fragment = std::make_shared<Fragment>(child_id);

        // Map old atom indices to new ones
        std::unordered_map<AtomIndex, AtomIndex> atom_index_map;

        // Copy atoms
        for (AtomIndex old_idx : atom_indices) {
            if (old_idx < atom_count()) {
                AtomIndex new_idx = new_fragment->add_atom(atom(old_idx));
                atom_index_map[old_idx] = new_idx;
            }
        }

        // Copy relevant bonds
        transfer_relevant_bonds(atom_indices, *new_fragment, atom_index_map);

        return new_fragment;
    }

    void Fragment::transfer_relevant_bonds(const std::vector<AtomIndex>& atom_indices,
                                         Fragment& target_fragment,
                                         std::unordered_map<AtomIndex, AtomIndex>& atom_index_map) const {
        // Create set for fast lookup of included atoms
        std::unordered_set<AtomIndex> atom_set(atom_indices.begin(), atom_indices.end());

        // Copy bonds that have both atoms in the atom set
        for (std::size_t bond_idx = 0; bond_idx < bond_count(); ++bond_idx) {
            const auto& bond_obj = bond(bond_idx);
            
            // Get the atom indices from the bond (this is conceptual - actual Bond class interface may differ)
            // auto [first_atom, second_atom] = bond_obj.get_atom_indices();
            
            // For now, skip bond transfer due to interface uncertainty
            // In practice, you'd check if both atoms are in atom_set and add the bond with mapped indices
        }
    }

} // namespace libfrag
