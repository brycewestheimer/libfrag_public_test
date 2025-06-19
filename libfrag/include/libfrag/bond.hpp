#pragma once
#ifndef LIBFRAG_BOND_HPP
#define LIBFRAG_BOND_HPP

#include <memory>
#include <optional>
#include <unordered_map>
#include <string>
#include <stdexcept>

namespace libfrag {

    // Forward declaration
    class Atom;

    /**
     * @brief Chemical bond class representing connection between two atoms
     * 
     * The Bond class represents a chemical bond between two atoms in a molecular system.
     * It stores references to the bonded atoms, bond properties such as order and length,
     * and additional quantum mechanical properties. The class is designed for compatibility
     * with quantum chemistry software and MolSSI tools.
     * 
     * @note Bonds maintain weak references to atoms to avoid circular dependencies.
     *       The atoms must remain valid throughout the bond's lifetime.
     */
    class Bond {
    public:
        /**
         * @brief Enumeration of common bond types
         */
        enum class BondType {
            SINGLE = 1,      ///< Single bond (bond order = 1)
            DOUBLE = 2,      ///< Double bond (bond order = 2) 
            TRIPLE = 3,      ///< Triple bond (bond order = 3)
            AROMATIC,        ///< Aromatic bond (delocalized)
            COORDINATE,      ///< Coordinate/dative bond
            HYDROGEN,        ///< Hydrogen bond (weak interaction)
            VAN_DER_WAALS,   ///< Van der Waals interaction
            IONIC,           ///< Ionic bond
            METALLIC,        ///< Metallic bond
            CUSTOM           ///< User-defined bond type
        };

        // Constructors
        Bond() = delete;  // Bonds must connect two atoms
        
        /**
         * @brief Construct bond between two atoms with specified type
         * 
         * @param first_atom Reference to first atom in the bond
         * @param second_atom Reference to second atom in the bond  
         * @param bond_type Type of chemical bond
         * 
         * @throws std::invalid_argument if atoms are the same object
         */
        Bond(const Atom& first_atom, const Atom& second_atom, BondType bond_type = BondType::SINGLE);
        
        /**
         * @brief Construct bond with explicit bond order
         * 
         * @param first_atom Reference to first atom in the bond
         * @param second_atom Reference to second atom in the bond
         * @param bond_order Fractional bond order (e.g., 1.5 for aromatic)
         */
        Bond(const Atom& first_atom, const Atom& second_atom, double bond_order);

        // Copy and move constructors
        Bond(const Bond& other) = default;
        Bond(Bond&& other) = default;
        Bond& operator=(const Bond& other) = default;
        Bond& operator=(Bond&& other) = default;
        
        // Destructor
        ~Bond() = default;

        // Accessors for atoms
        /**
         * @brief Get reference to first atom in bond
         * @return Const reference to first atom
         * @throws std::runtime_error if atom reference is invalid
         */
        const Atom& first_atom() const;
        
        /**
         * @brief Get reference to second atom in bond  
         * @return Const reference to second atom
         * @throws std::runtime_error if atom reference is invalid
         */
        const Atom& second_atom() const;
        
        /**
         * @brief Get the other atom in the bond given one atom
         * @param reference_atom One of the atoms in the bond
         * @return Reference to the other atom
         * @throws std::invalid_argument if reference_atom is not part of this bond
         */
        const Atom& get_other_atom(const Atom& reference_atom) const;
        
        /**
         * @brief Check if bond contains specified atom
         * @param atom_to_check Atom to search for
         * @return True if atom is part of this bond
         */
        bool contains_atom(const Atom& atom_to_check) const;

        // Bond properties
        BondType bond_type() const { return bond_type_; }
        double bond_order() const { return bond_order_; }
        
        /**
         * @brief Calculate current bond length from atomic coordinates
         * @return Bond length in Angstroms
         * @throws std::runtime_error if atom references are invalid
         */
        double bond_length() const;
        
        /**
         * @brief Get stored equilibrium bond length
         * @return Equilibrium bond length in Angstroms, or nullopt if not set
         */
        std::optional<double> equilibrium_bond_length() const { return equilibrium_bond_length_; }
        
        /**
         * @brief Get bond strength/energy
         * @return Bond dissociation energy in kcal/mol, or nullopt if not set
         */
        std::optional<double> bond_strength() const { return bond_strength_; }

        // Quantum properties
        /**
         * @brief Get all stored quantum properties
         * @return Map of property names to values
         */
        const std::unordered_map<std::string, double>& properties() const { return quantum_properties_; }
        
        /**
         * @brief Get specific quantum property value
         * @param property_name Name of the property to retrieve
         * @return Property value if found, nullopt otherwise
         */
        std::optional<double> get_property(const std::string& property_name) const;

        // Mutators
        void set_bond_type(BondType new_bond_type) { bond_type_ = new_bond_type; }
        void set_bond_order(double new_bond_order);
        void set_equilibrium_bond_length(double length_angstroms) { equilibrium_bond_length_ = length_angstroms; }
        void set_bond_strength(double energy_kcal_per_mol) { bond_strength_ = energy_kcal_per_mol; }
        
        /**
         * @brief Set quantum mechanical property
         * @param property_name Name of the property (e.g., "wiberg_bond_order", "mayer_bond_order")
         * @param property_value Numerical value of the property
         */
        void set_property(const std::string& property_name, double property_value);

        // Utility methods
        /**
         * @brief Check if bond is considered strong (covalent/ionic)
         * @return True for single, double, triple, or ionic bonds
         */
        bool is_strong_bond() const;
        
        /**
         * @brief Check if bond is considered weak (hydrogen, vdW)
         * @return True for hydrogen bonds and van der Waals interactions
         */
        bool is_weak_bond() const;
        
        /**
         * @brief Check if bond involves specified atom type
         * @param atomic_number Atomic number to check for
         * @return True if either atom has the specified atomic number
         */
        bool involves_atom_type(int atomic_number) const;
        
        /**
         * @brief Check if bond is polar (different electronegativity)
         * @return True if atoms have significantly different electronegativity
         */
        bool is_polar_bond() const;

        // String representations
        std::string to_string() const;
        std::string to_detailed_string() const;

        // Static utility methods
        /**
         * @brief Convert bond type enum to string
         * @param bond_type Bond type to convert
         * @return String representation of bond type
         */
        static std::string bond_type_to_string(BondType bond_type);
        
        /**
         * @brief Convert string to bond type enum
         * @param type_string String representation of bond type
         * @return Bond type enum, or nullopt if not recognized
         */
        static std::optional<BondType> string_to_bond_type(const std::string& type_string);
        
        /**
         * @brief Get typical bond order for bond type
         * @param bond_type Type of bond
         * @return Typical bond order (1.0 for single, 2.0 for double, etc.)
         */
        static double get_typical_bond_order(BondType bond_type);

        // Operators
        bool operator==(const Bond& other) const;
        bool operator!=(const Bond& other) const { return !(*this == other); }

    private:
        // Store weak references to atoms to avoid circular dependencies
        std::reference_wrapper<const Atom> first_atom_ref_;
        std::reference_wrapper<const Atom> second_atom_ref_;
        
        // Bond properties
        BondType bond_type_;
        double bond_order_;
        std::optional<double> equilibrium_bond_length_;  // Angstroms
        std::optional<double> bond_strength_;  // kcal/mol
        
        // Quantum mechanical properties (e.g., Wiberg bond orders, Mayer bond orders)
        std::unordered_map<std::string, double> quantum_properties_;
        
        // Helper methods
        void validate_bond_order(double order) const;
        double calculate_electronegativity_difference() const;
        
        // Static data for electronegativity values (Pauling scale)
        static const std::unordered_map<int, double> electronegativity_values_;
    };

    // Implementation of inline methods

    inline Bond::Bond(const Atom& first_atom, const Atom& second_atom, BondType bond_type)
        : first_atom_ref_(first_atom), second_atom_ref_(second_atom), 
          bond_type_(bond_type), bond_order_(get_typical_bond_order(bond_type)) {
        
        if (&first_atom == &second_atom) {
            throw std::invalid_argument("Cannot create bond between atom and itself");
        }
    }

    inline Bond::Bond(const Atom& first_atom, const Atom& second_atom, double bond_order)
        : first_atom_ref_(first_atom), second_atom_ref_(second_atom), 
          bond_type_(BondType::CUSTOM), bond_order_(bond_order) {
        
        if (&first_atom == &second_atom) {
            throw std::invalid_argument("Cannot create bond between atom and itself");
        }
        validate_bond_order(bond_order);
    }

    inline const Atom& Bond::first_atom() const {
        return first_atom_ref_.get();
    }

    inline const Atom& Bond::second_atom() const {
        return second_atom_ref_.get();
    }

    inline const Atom& Bond::get_other_atom(const Atom& reference_atom) const {
        if (&reference_atom == &first_atom()) {
            return second_atom();
        } else if (&reference_atom == &second_atom()) {
            return first_atom();
        } else {
            throw std::invalid_argument("Reference atom is not part of this bond");
        }
    }

    inline bool Bond::contains_atom(const Atom& atom_to_check) const {
        return (&atom_to_check == &first_atom()) || (&atom_to_check == &second_atom());
    }

    inline std::optional<double> Bond::get_property(const std::string& property_name) const {
        auto iterator = quantum_properties_.find(property_name);
        return (iterator != quantum_properties_.end()) ? 
               std::optional<double>(iterator->second) : std::nullopt;
    }

    inline void Bond::set_bond_order(double new_bond_order) {
        validate_bond_order(new_bond_order);
        bond_order_ = new_bond_order;
    }

    inline void Bond::set_property(const std::string& property_name, double property_value) {
        quantum_properties_[property_name] = property_value;
    }

    inline bool Bond::is_strong_bond() const {
        return bond_type_ == BondType::SINGLE || bond_type_ == BondType::DOUBLE || 
               bond_type_ == BondType::TRIPLE || bond_type_ == BondType::IONIC ||
               bond_type_ == BondType::COORDINATE;
    }

    inline bool Bond::is_weak_bond() const {
        return bond_type_ == BondType::HYDROGEN || bond_type_ == BondType::VAN_DER_WAALS;
    }

    inline void Bond::validate_bond_order(double order) const {
        if (order < 0.0 || order > 4.0) {
            throw std::invalid_argument("Bond order must be between 0.0 and 4.0");
        }
    }

    // Static data for electronegativity values (Pauling scale)
    inline const std::unordered_map<int, double> Bond::electronegativity_values_ = {
        {1, 2.20}, {2, 0.0}, {3, 0.98}, {4, 1.57}, {5, 2.04}, {6, 2.55}, {7, 3.04}, {8, 3.44},
        {9, 3.98}, {10, 0.0}, {11, 0.93}, {12, 1.31}, {13, 1.61}, {14, 1.90}, {15, 2.19}, {16, 2.58},
        {17, 3.16}, {18, 0.0}, {19, 0.82}, {20, 1.00}, {21, 1.36}, {22, 1.54}, {23, 1.63}, {24, 1.66},
        {25, 1.55}, {26, 1.83}, {27, 1.88}, {28, 1.91}, {29, 1.90}, {30, 1.65}, {31, 1.81}, {32, 2.01},
        {33, 2.18}, {34, 2.55}, {35, 2.96}, {36, 3.00}
    };

} // namespace libfrag

#endif // LIBFRAG_BOND_HPP
