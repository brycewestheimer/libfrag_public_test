#include "libfrag/molecule.hpp"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <set>
#include <queue>

namespace libfrag {

    Molecule::Molecule(const std::vector<Atom>& atom_list) {
        atoms_.reserve(atom_list.size());
        for (const auto& atom : atom_list) {
            atoms_.emplace_back(std::make_unique<Atom>(atom));
        }
    }

    Molecule::Molecule(const std::vector<Atom>& atom_list,
                      const std::vector<std::pair<AtomIndex, AtomIndex>>& bond_connectivity,
                      const std::vector<Bond::BondType>& bond_types) : Molecule(atom_list) {
        
        // Add bonds based on connectivity
        for (std::size_t bond_index = 0; bond_index < bond_connectivity.size(); ++bond_index) {
            auto [first_atom_index, second_atom_index] = bond_connectivity[bond_index];
            
            Bond::BondType bond_type = Bond::BondType::SINGLE;  // Default
            if (bond_index < bond_types.size()) {
                bond_type = bond_types[bond_index];
            }
            
            add_bond(first_atom_index, second_atom_index, bond_type);
        }
    }

    Molecule::Molecule(const Molecule& other) 
        : molecular_multiplicity_(other.molecular_multiplicity_),
          molecular_properties_(other.molecular_properties_) {
        
        // Deep copy atoms
        atoms_.reserve(other.atoms_.size());
        for (const auto& atom_ptr : other.atoms_) {
            atoms_.emplace_back(std::make_unique<Atom>(*atom_ptr));
        }
        
        // Deep copy bonds with new atom references
        bonds_.reserve(other.bonds_.size());
        for (const auto& bond_ptr : other.bonds_) {
            const Bond& original_bond = *bond_ptr;
            
            // Find indices of the bonded atoms in original molecule
            AtomIndex first_atom_index = 0, second_atom_index = 0;
            bool first_found = false, second_found = false;
            
            for (std::size_t atom_index = 0; atom_index < other.atoms_.size(); ++atom_index) {
                if (&original_bond.first_atom() == other.atoms_[atom_index].get()) {
                    first_atom_index = atom_index;
                    first_found = true;
                }
                if (&original_bond.second_atom() == other.atoms_[atom_index].get()) {
                    second_atom_index = atom_index;
                    second_found = true;
                }
                if (first_found && second_found) break;
            }
            
            // Create new bond with copied atoms
            auto new_bond = std::make_unique<Bond>(*atoms_[first_atom_index], 
                                                  *atoms_[second_atom_index], 
                                                  original_bond.bond_order());
            new_bond->set_bond_type(original_bond.bond_type());
            
            // Copy bond properties
            for (const auto& property_pair : original_bond.properties()) {
                new_bond->set_property(property_pair.first, property_pair.second);
            }
            
            bonds_.emplace_back(std::move(new_bond));
        }
    }

    Molecule& Molecule::operator=(const Molecule& other) {
        if (this != &other) {
            Molecule temp_copy(other);  // Create temporary copy
            *this = std::move(temp_copy);  // Move assign from temporary
        }
        return *this;
    }

    AtomIndex Molecule::add_atom(const Atom& new_atom) {
        atoms_.emplace_back(std::make_unique<Atom>(new_atom));
        return atoms_.size() - 1;
    }

    void Molecule::remove_atom(AtomIndex atom_index) {
        validate_atom_index(atom_index);
        
        // Remove all bonds involving this atom
        remove_bonds_involving_atom(atom_index);
        
        // Remove the atom itself
        atoms_.erase(atoms_.begin() + atom_index);
        
        // Update bond references to account for shifted indices
        // This is complex - bonds store references, so we need to rebuild them
        // For now, we'll keep it simple and require users to be careful with indices
    }

    BondIndex Molecule::add_bond(AtomIndex first_atom_index, AtomIndex second_atom_index, 
                                Bond::BondType bond_type) {
        validate_atom_index(first_atom_index);
        validate_atom_index(second_atom_index);
        
        if (bond_exists(first_atom_index, second_atom_index)) {
            throw std::invalid_argument("Bond already exists between atoms " + 
                                      std::to_string(first_atom_index) + " and " + 
                                      std::to_string(second_atom_index));
        }
        
        bonds_.emplace_back(std::make_unique<Bond>(*atoms_[first_atom_index], 
                                                  *atoms_[second_atom_index], 
                                                  bond_type));
        return bonds_.size() - 1;
    }

    BondIndex Molecule::add_bond(AtomIndex first_atom_index, AtomIndex second_atom_index, 
                                double bond_order) {
        validate_atom_index(first_atom_index);
        validate_atom_index(second_atom_index);
        
        if (bond_exists(first_atom_index, second_atom_index)) {
            throw std::invalid_argument("Bond already exists between atoms " + 
                                      std::to_string(first_atom_index) + " and " + 
                                      std::to_string(second_atom_index));
        }
        
        bonds_.emplace_back(std::make_unique<Bond>(*atoms_[first_atom_index], 
                                                  *atoms_[second_atom_index], 
                                                  bond_order));
        return bonds_.size() - 1;
    }

    void Molecule::remove_bond(BondIndex bond_index) {
        validate_bond_index(bond_index);
        bonds_.erase(bonds_.begin() + bond_index);
    }

    std::optional<BondIndex> Molecule::find_bond(AtomIndex first_atom_index, 
                                                AtomIndex second_atom_index) const {
        validate_atom_index(first_atom_index);
        validate_atom_index(second_atom_index);
        
        for (std::size_t bond_index = 0; bond_index < bonds_.size(); ++bond_index) {
            const Bond& current_bond = *bonds_[bond_index];
            if (current_bond.contains_atom(*atoms_[first_atom_index]) &&
                current_bond.contains_atom(*atoms_[second_atom_index])) {
                return bond_index;
            }
        }
        return std::nullopt;
    }

    std::vector<BondIndex> Molecule::bonds_for_atom(AtomIndex atom_index) const {
        validate_atom_index(atom_index);
        
        std::vector<BondIndex> result_bonds;
        const Atom& target_atom = *atoms_[atom_index];
        
        for (std::size_t bond_index = 0; bond_index < bonds_.size(); ++bond_index) {
            if (bonds_[bond_index]->contains_atom(target_atom)) {
                result_bonds.push_back(bond_index);
            }
        }
        
        return result_bonds;
    }

    double Molecule::molecular_mass() const {
        double total_mass = 0.0;
        for (const auto& atom_ptr : atoms_) {
            total_mass += atom_ptr->atomic_mass();
        }
        return total_mass;
    }

    std::array<double, 3> Molecule::center_of_mass() const {
        if (atoms_.empty()) {
            return {0.0, 0.0, 0.0};
        }
        
        double total_mass = 0.0;
        std::array<double, 3> weighted_position = {0.0, 0.0, 0.0};
        
        for (const auto& atom_ptr : atoms_) {
            double atom_mass = atom_ptr->atomic_mass();
            total_mass += atom_mass;
            weighted_position[0] += atom_mass * atom_ptr->x();
            weighted_position[1] += atom_mass * atom_ptr->y();
            weighted_position[2] += atom_mass * atom_ptr->z();
        }
        
        weighted_position[0] /= total_mass;
        weighted_position[1] /= total_mass;
        weighted_position[2] /= total_mass;
        
        return weighted_position;
    }

    double Molecule::molecular_charge() const {
        double total_charge = 0.0;
        for (const auto& atom_ptr : atoms_) {
            total_charge += atom_ptr->formal_charge();
        }
        return total_charge;
    }

    void Molecule::translate(const std::array<double, 3>& displacement_vector) {
        for (auto& atom_ptr : atoms_) {
            auto current_coords = atom_ptr->coordinates();
            atom_ptr->set_coordinates(current_coords[0] + displacement_vector[0],
                                    current_coords[1] + displacement_vector[1],
                                    current_coords[2] + displacement_vector[2]);
        }
    }

    void Molecule::center_at_origin() {
        auto com = center_of_mass();
        translate({-com[0], -com[1], -com[2]});
    }

    std::array<double, 6> Molecule::bounding_box() const {
        if (atoms_.empty()) {
            return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        }
        
        double min_x = atoms_[0]->x(), max_x = atoms_[0]->x();
        double min_y = atoms_[0]->y(), max_y = atoms_[0]->y();
        double min_z = atoms_[0]->z(), max_z = atoms_[0]->z();
        
        for (const auto& atom_ptr : atoms_) {
            min_x = std::min(min_x, atom_ptr->x());
            max_x = std::max(max_x, atom_ptr->x());
            min_y = std::min(min_y, atom_ptr->y());
            max_y = std::max(max_y, atom_ptr->y());
            min_z = std::min(min_z, atom_ptr->z());
            max_z = std::max(max_z, atom_ptr->z());
        }
        
        return {min_x, min_y, min_z, max_x, max_y, max_z};
    }

    double Molecule::molecular_diameter() const {
        if (atoms_.size() < 2) {
            return 0.0;
        }
        
        double max_distance = 0.0;
        for (std::size_t i = 0; i < atoms_.size(); ++i) {
            for (std::size_t j = i + 1; j < atoms_.size(); ++j) {
                double distance = atoms_[i]->distance_to(*atoms_[j]);
                max_distance = std::max(max_distance, distance);
            }
        }
        
        return max_distance;
    }

    std::size_t Molecule::coordination_number(AtomIndex atom_index) const {
        return bonds_for_atom(atom_index).size();
    }

    std::vector<AtomIndex> Molecule::bonded_atoms(AtomIndex atom_index) const {
        validate_atom_index(atom_index);
        
        std::vector<AtomIndex> bonded_atom_indices;
        const Atom& central_atom = *atoms_[atom_index];
        
        for (const auto& bond_ptr : bonds_) {
            if (bond_ptr->contains_atom(central_atom)) {
                // Find the other atom in the bond
                const Atom& other_atom = bond_ptr->get_other_atom(central_atom);
                
                // Find its index
                for (std::size_t other_index = 0; other_index < atoms_.size(); ++other_index) {
                    if (atoms_[other_index].get() == &other_atom) {
                        bonded_atom_indices.push_back(other_index);
                        break;
                    }
                }
            }
        }
        
        return bonded_atom_indices;
    }

    bool Molecule::is_connected() const {
        if (atoms_.size() <= 1) {
            return true;  // Single atom or empty molecule is considered connected
        }
        
        std::vector<bool> visited_atoms(atoms_.size(), false);
        std::vector<AtomIndex> current_fragment;
        
        // Start DFS from first atom
        depth_first_search_connectivity(0, visited_atoms, current_fragment);
        
        // Check if all atoms were visited
        return current_fragment.size() == atoms_.size();
    }

    std::vector<Molecule> Molecule::get_fragments() const {
        std::vector<Molecule> molecular_fragments;
        std::vector<bool> visited_atoms(atoms_.size(), false);
        
        for (std::size_t start_atom_index = 0; start_atom_index < atoms_.size(); ++start_atom_index) {
            if (!visited_atoms[start_atom_index]) {
                std::vector<AtomIndex> fragment_atom_indices;
                depth_first_search_connectivity(start_atom_index, visited_atoms, fragment_atom_indices);
                
                // Create molecule for this fragment
                std::vector<Atom> fragment_atoms;
                for (AtomIndex atom_index : fragment_atom_indices) {
                    fragment_atoms.push_back(*atoms_[atom_index]);
                }
                
                // Find bonds within this fragment
                std::vector<std::pair<AtomIndex, AtomIndex>> fragment_bonds;
                std::vector<Bond::BondType> fragment_bond_types;
                
                for (const auto& bond_ptr : bonds_) {
                    // Check if both atoms of this bond are in the current fragment
                    bool first_atom_in_fragment = false, second_atom_in_fragment = false;
                    AtomIndex first_local_index = 0, second_local_index = 0;
                    
                    for (std::size_t local_index = 0; local_index < fragment_atom_indices.size(); ++local_index) {
                        AtomIndex global_index = fragment_atom_indices[local_index];
                        if (atoms_[global_index].get() == &bond_ptr->first_atom()) {
                            first_atom_in_fragment = true;
                            first_local_index = local_index;
                        }
                        if (atoms_[global_index].get() == &bond_ptr->second_atom()) {
                            second_atom_in_fragment = true;
                            second_local_index = local_index;
                        }
                    }
                    
                    if (first_atom_in_fragment && second_atom_in_fragment) {
                        fragment_bonds.emplace_back(first_local_index, second_local_index);
                        fragment_bond_types.push_back(bond_ptr->bond_type());
                    }
                }
                
                molecular_fragments.emplace_back(fragment_atoms, fragment_bonds, fragment_bond_types);
            }
        }
        
        return molecular_fragments;
    }

    std::vector<int> Molecule::atomic_numbers() const {
        std::vector<int> atomic_number_list;
        atomic_number_list.reserve(atoms_.size());
        
        for (const auto& atom_ptr : atoms_) {
            atomic_number_list.push_back(atom_ptr->atomic_number());
        }
        
        return atomic_number_list;
    }

    std::vector<std::string> Molecule::element_symbols() const {
        std::vector<std::string> symbol_list;
        symbol_list.reserve(atoms_.size());
        
        for (const auto& atom_ptr : atoms_) {
            symbol_list.push_back(atom_ptr->symbol());
        }
        
        return symbol_list;
    }

    xt::xarray<double> Molecule::geometry_array() const {
        xt::xarray<double> geometry_matrix = xt::zeros<double>({atoms_.size(), 3});
        
        for (std::size_t atom_index = 0; atom_index < atoms_.size(); ++atom_index) {
            const auto& atom_coords = atoms_[atom_index]->coordinates();
            geometry_matrix(atom_index, 0) = atom_coords[0];
            geometry_matrix(atom_index, 1) = atom_coords[1];
            geometry_matrix(atom_index, 2) = atom_coords[2];
        }
        
        return xt::flatten(geometry_matrix);
    }

    std::vector<std::array<double, 3>> Molecule::geometry_vector() const {
        std::vector<std::array<double, 3>> coordinate_list;
        coordinate_list.reserve(atoms_.size());
        
        for (const auto& atom_ptr : atoms_) {
            coordinate_list.push_back(atom_ptr->coordinates());
        }
        
        return coordinate_list;
    }

    xt::xarray<double> Molecule::connectivity_matrix() const {
        std::size_t atom_count = atoms_.size();
        xt::xarray<double> connectivity = xt::zeros<double>({atom_count, atom_count});
        
        for (const auto& bond_ptr : bonds_) {
            // Find indices of bonded atoms
            AtomIndex first_index = 0, second_index = 0;
            bool first_found = false, second_found = false;
            
            for (std::size_t atom_index = 0; atom_index < atoms_.size(); ++atom_index) {
                if (atoms_[atom_index].get() == &bond_ptr->first_atom()) {
                    first_index = atom_index;
                    first_found = true;
                }
                if (atoms_[atom_index].get() == &bond_ptr->second_atom()) {
                    second_index = atom_index;
                    second_found = true;
                }
                if (first_found && second_found) break;
            }
            
            if (first_found && second_found) {
                connectivity(first_index, second_index) = bond_ptr->bond_order();
                connectivity(second_index, first_index) = bond_ptr->bond_order();
            }
        }
        
        return connectivity;
    }

    std::string Molecule::to_xyz_string(const std::string& comment_line) const {
        std::stringstream xyz_stream;
        xyz_stream << atoms_.size() << "\n";
        xyz_stream << comment_line << "\n";
        
        for (const auto& atom_ptr : atoms_) {
            xyz_stream << atom_ptr->to_xyz_string() << "\n";
        }
        
        return xyz_stream.str();
    }

    std::string Molecule::to_detailed_string() const {
        std::stringstream detail_stream;
        detail_stream << std::fixed << std::setprecision(3);
        
        detail_stream << "Molecule Information:\n";
        detail_stream << "  Atoms: " << atoms_.size() << "\n";
        detail_stream << "  Bonds: " << bonds_.size() << "\n";
        detail_stream << "  Molecular mass: " << molecular_mass() << " amu\n";
        detail_stream << "  Molecular charge: " << molecular_charge() << "\n";
        detail_stream << "  Multiplicity: " << molecular_multiplicity_ << "\n";
        
        auto com = center_of_mass();
        detail_stream << "  Center of mass: (" << com[0] << ", " << com[1] << ", " << com[2] << ")\n";
        
        detail_stream << "\nAtoms:\n";
        for (std::size_t atom_index = 0; atom_index < atoms_.size(); ++atom_index) {
            detail_stream << "  " << atom_index << ": " << atoms_[atom_index]->to_xyz_string() << "\n";
        }
        
        if (!bonds_.empty()) {
            detail_stream << "\nBonds:\n";
            for (std::size_t bond_index = 0; bond_index < bonds_.size(); ++bond_index) {
                detail_stream << "  " << bond_index << ": " << bonds_[bond_index]->to_string() << "\n";
            }
        }
        
        return detail_stream.str();
    }

    std::string Molecule::to_string() const {
        if (atoms_.empty()) {
            return "Empty molecule";
        }
        
        // Count atoms by element
        std::map<std::string, int> element_counts;
        for (const auto& atom_ptr : atoms_) {
            element_counts[atom_ptr->symbol()]++;
        }
        
        // Build molecular formula
        std::stringstream formula_stream;
        for (const auto& element_pair : element_counts) {
            formula_stream << element_pair.first;
            if (element_pair.second > 1) {
                formula_stream << element_pair.second;
            }
        }
        
        return formula_stream.str() + " (" + std::to_string(atoms_.size()) + " atoms, " + 
               std::to_string(bonds_.size()) + " bonds)";
    }

    bool Molecule::operator==(const Molecule& other) const {
        // Two molecules are equal if they have the same atoms and bonds
        if (atoms_.size() != other.atoms_.size() || bonds_.size() != other.bonds_.size()) {
            return false;
        }
        
        // This is a simplified equality check - in practice, we might need
        // more sophisticated graph isomorphism checking
        for (std::size_t atom_index = 0; atom_index < atoms_.size(); ++atom_index) {
            if (*atoms_[atom_index] != *other.atoms_[atom_index]) {
                return false;
            }
        }
        
        return true;
    }

    bool Molecule::bond_exists(AtomIndex first_index, AtomIndex second_index) const {
        return find_bond(first_index, second_index).has_value();
    }

    void Molecule::remove_bonds_involving_atom(AtomIndex atom_index) {
        const Atom& target_atom = *atoms_[atom_index];
        
        bonds_.erase(
            std::remove_if(bonds_.begin(), bonds_.end(),
                          [&target_atom](const std::unique_ptr<Bond>& bond_ptr) {
                              return bond_ptr->contains_atom(target_atom);
                          }),
            bonds_.end()
        );
    }

    void Molecule::depth_first_search_connectivity(AtomIndex current_atom, 
                                                  std::vector<bool>& visited_atoms,
                                                  std::vector<AtomIndex>& current_fragment) const {
        visited_atoms[current_atom] = true;
        current_fragment.push_back(current_atom);
        
        // Visit all bonded atoms
        auto bonded_atom_indices = bonded_atoms(current_atom);
        for (AtomIndex bonded_atom : bonded_atom_indices) {
            if (!visited_atoms[bonded_atom]) {
                depth_first_search_connectivity(bonded_atom, visited_atoms, current_fragment);
            }
        }
    }

    // Static methods
    Molecule Molecule::from_xyz_string(const std::string& xyz_content) {
        std::istringstream input_stream(xyz_content);
        std::string line_buffer;
        
        // Read number of atoms
        std::getline(input_stream, line_buffer);
        int atom_count = std::stoi(line_buffer);
        
        // Skip comment line
        std::getline(input_stream, line_buffer);
        
        // Read atoms
        std::vector<Atom> atom_list;
        atom_list.reserve(atom_count);
        
        for (int atom_index = 0; atom_index < atom_count; ++atom_index) {
            std::getline(input_stream, line_buffer);
            std::istringstream line_stream(line_buffer);
            
            std::string element_symbol;
            double x_coord, y_coord, z_coord;
            
            line_stream >> element_symbol >> x_coord >> y_coord >> z_coord;
            atom_list.emplace_back(element_symbol, x_coord, y_coord, z_coord);
        }
        
        return Molecule(atom_list);
    }

    Molecule Molecule::create_linear_molecule(const std::vector<std::string>& element_symbols, 
                                            double bond_length) {
        if (element_symbols.empty()) {
            return Molecule();
        }
        
        std::vector<Atom> linear_atoms;
        linear_atoms.reserve(element_symbols.size());
        
        for (std::size_t atom_index = 0; atom_index < element_symbols.size(); ++atom_index) {
            double x_position = static_cast<double>(atom_index) * bond_length;
            linear_atoms.emplace_back(element_symbols[atom_index], x_position, 0.0, 0.0);
        }
        
        // Create bonds between adjacent atoms
        std::vector<std::pair<AtomIndex, AtomIndex>> bond_connectivity;
        for (std::size_t atom_index = 0; atom_index < element_symbols.size() - 1; ++atom_index) {
            bond_connectivity.emplace_back(atom_index, atom_index + 1);
        }
        
        return Molecule(linear_atoms, bond_connectivity);
    }

} // namespace libfrag
