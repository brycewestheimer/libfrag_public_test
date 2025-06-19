#include "libfrag/bond.hpp"
#include "libfrag/atom.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace libfrag {

    double Bond::bond_length() const {
        try {
            return first_atom().distance_to(second_atom());
        } catch (const std::exception& e) {
            throw std::runtime_error("Cannot calculate bond length: " + std::string(e.what()));
        }
    }

    bool Bond::involves_atom_type(int atomic_number) const {
        return (first_atom().atomic_number() == atomic_number) || 
               (second_atom().atomic_number() == atomic_number);
    }

    bool Bond::is_polar_bond() const {
        const double electronegativity_threshold = 0.4;  // Pauling scale difference
        double electronegativity_difference = calculate_electronegativity_difference();
        return electronegativity_difference > electronegativity_threshold;
    }

    std::string Bond::to_string() const {
        std::stringstream stream_buffer;
        stream_buffer << first_atom().symbol() << "-" << second_atom().symbol() 
                     << " (" << bond_type_to_string(bond_type_) 
                     << ", order=" << std::fixed << std::setprecision(2) << bond_order_ << ")";
        return stream_buffer.str();
    }

    std::string Bond::to_detailed_string() const {
        std::stringstream stream_buffer;
        stream_buffer << std::fixed << std::setprecision(3);
        
        stream_buffer << "Bond: " << first_atom().symbol() << "-" << second_atom().symbol() << "\n";
        stream_buffer << "  Type: " << bond_type_to_string(bond_type_) << "\n";
        stream_buffer << "  Order: " << bond_order_ << "\n";
        stream_buffer << "  Length: " << bond_length() << " Angstrom\n";
        
        if (equilibrium_bond_length_) {
            stream_buffer << "  Equilibrium Length: " << *equilibrium_bond_length_ << " Angstrom\n";
        }
        
        if (bond_strength_) {
            stream_buffer << "  Bond Strength: " << *bond_strength_ << " kcal/mol\n";
        }
        
        if (!quantum_properties_.empty()) {
            stream_buffer << "  Quantum Properties:\n";
            for (const auto& property_pair : quantum_properties_) {
                stream_buffer << "    " << property_pair.first << ": " << property_pair.second << "\n";
            }
        }
        
        return stream_buffer.str();
    }

    std::string Bond::bond_type_to_string(BondType bond_type) {
        switch (bond_type) {
            case BondType::SINGLE:       return "Single";
            case BondType::DOUBLE:       return "Double";
            case BondType::TRIPLE:       return "Triple";
            case BondType::AROMATIC:     return "Aromatic";
            case BondType::COORDINATE:   return "Coordinate";
            case BondType::HYDROGEN:     return "Hydrogen";
            case BondType::VAN_DER_WAALS: return "van der Waals";
            case BondType::IONIC:        return "Ionic";
            case BondType::METALLIC:     return "Metallic";
            case BondType::CUSTOM:       return "Custom";
            default:                     return "Unknown";
        }
    }

    std::optional<Bond::BondType> Bond::string_to_bond_type(const std::string& type_string) {
        // Convert to lowercase for case-insensitive comparison
        std::string lowercase_string = type_string;
        std::transform(lowercase_string.begin(), lowercase_string.end(), 
                      lowercase_string.begin(), ::tolower);
        
        if (lowercase_string == "single" || lowercase_string == "1") {
            return BondType::SINGLE;
        } else if (lowercase_string == "double" || lowercase_string == "2") {
            return BondType::DOUBLE;
        } else if (lowercase_string == "triple" || lowercase_string == "3") {
            return BondType::TRIPLE;
        } else if (lowercase_string == "aromatic") {
            return BondType::AROMATIC;
        } else if (lowercase_string == "coordinate" || lowercase_string == "dative") {
            return BondType::COORDINATE;
        } else if (lowercase_string == "hydrogen" || lowercase_string == "h-bond") {
            return BondType::HYDROGEN;
        } else if (lowercase_string == "van der waals" || lowercase_string == "vdw") {
            return BondType::VAN_DER_WAALS;
        } else if (lowercase_string == "ionic") {
            return BondType::IONIC;
        } else if (lowercase_string == "metallic") {
            return BondType::METALLIC;
        } else if (lowercase_string == "custom") {
            return BondType::CUSTOM;
        } else {
            return std::nullopt;
        }
    }

    double Bond::get_typical_bond_order(BondType bond_type) {
        switch (bond_type) {
            case BondType::SINGLE:       return 1.0;
            case BondType::DOUBLE:       return 2.0;
            case BondType::TRIPLE:       return 3.0;
            case BondType::AROMATIC:     return 1.5;  // Typical aromatic bond order
            case BondType::COORDINATE:   return 1.0;  // Usually single bond strength
            case BondType::HYDROGEN:     return 0.1;  // Weak interaction
            case BondType::VAN_DER_WAALS: return 0.05; // Very weak interaction
            case BondType::IONIC:        return 1.0;  // Treated as single bond equivalent
            case BondType::METALLIC:     return 0.5;  // Variable, but typically less than single
            case BondType::CUSTOM:       return 1.0;  // Default value
            default:                     return 1.0;
        }
    }

    bool Bond::operator==(const Bond& other) const {
        // Two bonds are equal if they connect the same atoms (order doesn't matter)
        // and have the same bond type and order
        const double bond_order_tolerance = 1e-6;
        
        bool same_atoms = ((&first_atom() == &other.first_atom() && &second_atom() == &other.second_atom()) ||
                          (&first_atom() == &other.second_atom() && &second_atom() == &other.first_atom()));
        
        bool same_type_and_order = (bond_type_ == other.bond_type_) &&
                                  (std::abs(bond_order_ - other.bond_order_) < bond_order_tolerance);
        
        return same_atoms && same_type_and_order;
    }

    double Bond::calculate_electronegativity_difference() const {
        int first_atom_number = first_atom().atomic_number();
        int second_atom_number = second_atom().atomic_number();
        
        auto first_electronegativity_iterator = electronegativity_values_.find(first_atom_number);
        auto second_electronegativity_iterator = electronegativity_values_.find(second_atom_number);
        
        if (first_electronegativity_iterator == electronegativity_values_.end() ||
            second_electronegativity_iterator == electronegativity_values_.end()) {
            // Return 0 if electronegativity data not available
            return 0.0;
        }
        
        double first_electronegativity = first_electronegativity_iterator->second;
        double second_electronegativity = second_electronegativity_iterator->second;
        
        return std::abs(first_electronegativity - second_electronegativity);
    }

} // namespace libfrag
