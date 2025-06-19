#pragma once
#ifndef LIBFRAG_MOLECULE_HPP
#define LIBFRAG_MOLECULE_HPP

#include "libfrag/atom.hpp"
#include "libfrag/bond.hpp"
#include <vector>
#include <memory>
#include <unordered_map>
#include <string>
#include <optional>
#include <array>
#include <stdexcept>
#include "xtensor/xarray.hpp"

namespace libfrag {

    /**
     * @brief Comprehensive molecular system class for quantum chemistry
     * 
     * The Molecule class represents a complete molecular system containing atoms and bonds.
     * It provides methods for molecular property calculation, geometry manipulation,
     * and compatibility with quantum chemistry software packages. The class maintains
     * unique atom instances and manages bonds between them using references.
     * 
     * Key features:
     * - Unique atom management with automatic indexing
     * - Bond connectivity tracking with reference-based implementation
     * - Molecular property calculation (center of mass, moments of inertia, etc.)
     * - Quantum chemistry software compatibility (QCElemental, QCEngine)
     * - Molecular geometry operations (translation, rotation, optimization tracking)
     * - Fragment analysis and manipulation capabilities
     */
    class Molecule {
    public:
        // Type aliases for clarity
        using AtomContainer = std::vector<std::unique_ptr<Atom>>;
        using BondContainer = std::vector<std::unique_ptr<Bond>>;
        using AtomIndex = std::size_t;
        using BondIndex = std::size_t;

        // Constructors
        /**
         * @brief Default constructor creates empty molecule
         */
        Molecule() = default;
        
        /**
         * @brief Construct molecule from list of atoms (no bonds)
         * @param atom_list Vector of atoms to copy into molecule
         */
        explicit Molecule(const std::vector<Atom>& atom_list);
        
        /**
         * @brief Construct molecule with atoms and connectivity
         * @param atom_list Vector of atoms to copy
         * @param bond_connectivity Vector of atom index pairs representing bonds
         * @param bond_types Optional bond types (defaults to single bonds)
         */
        Molecule(const std::vector<Atom>& atom_list,
                const std::vector<std::pair<AtomIndex, AtomIndex>>& bond_connectivity,
                const std::vector<Bond::BondType>& bond_types = {});

        // Copy constructor and assignment operator
        Molecule(const Molecule& other);
        Molecule& operator=(const Molecule& other);
        
        // Move constructor and assignment operator  
        Molecule(Molecule&& other) = default;
        Molecule& operator=(Molecule&& other) = default;
        
        // Destructor
        ~Molecule() = default;

        // Basic properties
        /**
         * @brief Get number of atoms in molecule
         * @return Number of atoms
         */
        std::size_t atom_count() const { return atoms_.size(); }
        
        /**
         * @brief Get number of bonds in molecule
         * @return Number of bonds
         */
        std::size_t bond_count() const { return bonds_.size(); }
        
        /**
         * @brief Check if molecule is empty (no atoms)
         * @return True if molecule contains no atoms
         */
        bool is_empty() const { return atoms_.empty(); }

        // Atom access and manipulation
        /**
         * @brief Get atom by index
         * @param atom_index Zero-based index of atom
         * @return Const reference to atom
         * @throws std::out_of_range if index is invalid
         */
        const Atom& atom(AtomIndex atom_index) const;
        
        /**
         * @brief Get mutable atom by index
         * @param atom_index Zero-based index of atom
         * @return Reference to atom
         * @throws std::out_of_range if index is invalid
         */
        Atom& atom(AtomIndex atom_index);
        
        /**
         * @brief Add new atom to molecule
         * @param new_atom Atom to add (will be copied)
         * @return Index of newly added atom
         */
        AtomIndex add_atom(const Atom& new_atom);
        
        /**
         * @brief Add new atom constructed in-place
         * @param atomic_number Atomic number of new atom
         * @param x X coordinate in Angstroms
         * @param y Y coordinate in Angstroms  
         * @param z Z coordinate in Angstroms
         * @return Index of newly added atom
         */
        AtomIndex add_atom(int atomic_number, double x, double y, double z);
        
        /**
         * @brief Remove atom and all associated bonds
         * @param atom_index Index of atom to remove
         * @throws std::out_of_range if index is invalid
         */
        void remove_atom(AtomIndex atom_index);

        // Bond access and manipulation
        /**
         * @brief Get bond by index
         * @param bond_index Zero-based index of bond
         * @return Const reference to bond
         * @throws std::out_of_range if index is invalid
         */
        const Bond& bond(BondIndex bond_index) const;
        
        /**
         * @brief Get mutable bond by index
         * @param bond_index Zero-based index of bond
         * @return Reference to bond
         * @throws std::out_of_range if index is invalid
         */
        Bond& bond(BondIndex bond_index);
        
        /**
         * @brief Add bond between two atoms
         * @param first_atom_index Index of first atom
         * @param second_atom_index Index of second atom
         * @param bond_type Type of bond to create
         * @return Index of newly created bond
         * @throws std::out_of_range if atom indices are invalid
         * @throws std::invalid_argument if bond already exists
         */
        BondIndex add_bond(AtomIndex first_atom_index, AtomIndex second_atom_index, 
                          Bond::BondType bond_type = Bond::BondType::SINGLE);
        
        /**
         * @brief Add bond with explicit bond order
         * @param first_atom_index Index of first atom
         * @param second_atom_index Index of second atom
         * @param bond_order Fractional bond order
         * @return Index of newly created bond
         */
        BondIndex add_bond(AtomIndex first_atom_index, AtomIndex second_atom_index, double bond_order);
        
        /**
         * @brief Remove bond between atoms
         * @param bond_index Index of bond to remove
         * @throws std::out_of_range if index is invalid
         */
        void remove_bond(BondIndex bond_index);
        
        /**
         * @brief Find bond connecting two atoms
         * @param first_atom_index Index of first atom
         * @param second_atom_index Index of second atom
         * @return Bond index if found, nullopt otherwise
         */
        std::optional<BondIndex> find_bond(AtomIndex first_atom_index, AtomIndex second_atom_index) const;
        
        /**
         * @brief Get all bonds involving specified atom
         * @param atom_index Index of atom
         * @return Vector of bond indices
         */
        std::vector<BondIndex> bonds_for_atom(AtomIndex atom_index) const;

        // Molecular properties
        /**
         * @brief Calculate total molecular mass
         * @return Molecular mass in atomic mass units
         */
        double molecular_mass() const;
        
        /**
         * @brief Calculate center of mass
         * @return Center of mass coordinates in Angstroms
         */
        std::array<double, 3> center_of_mass() const;
        
        /**
         * @brief Calculate molecular charge (sum of atomic formal charges)
         * @return Total molecular charge
         */
        double molecular_charge() const;
        
        /**
         * @brief Get/set molecular multiplicity (2S + 1)
         * @return Spin multiplicity
         */
        int molecular_multiplicity() const { return molecular_multiplicity_; }
        void set_molecular_multiplicity(int multiplicity) { molecular_multiplicity_ = multiplicity; }

        // Geometry operations
        /**
         * @brief Translate entire molecule by vector
         * @param displacement_vector Translation vector [dx, dy, dz] in Angstroms
         */
        void translate(const std::array<double, 3>& displacement_vector);
        
        /**
         * @brief Translate molecule to place center of mass at origin
         */
        void center_at_origin();
        
        /**
         * @brief Get bounding box of molecular geometry
         * @return Array of [min_x, min_y, min_z, max_x, max_y, max_z]
         */
        std::array<double, 6> bounding_box() const;
        
        /**
         * @brief Calculate maximum distance between any two atoms
         * @return Maximum interatomic distance in Angstroms
         */
        double molecular_diameter() const;

        // Connectivity analysis
        /**
         * @brief Get coordination number (number of bonds) for atom
         * @param atom_index Index of atom
         * @return Number of bonds to this atom
         */
        std::size_t coordination_number(AtomIndex atom_index) const;
        
        /**
         * @brief Get atoms bonded to specified atom
         * @param atom_index Index of central atom
         * @return Vector of bonded atom indices
         */
        std::vector<AtomIndex> bonded_atoms(AtomIndex atom_index) const;
        
        /**
         * @brief Check if molecule is connected (single component)
         * @return True if all atoms are connected through bonds
         */
        bool is_connected() const;
        
        /**
         * @brief Get separate molecular fragments
         * @return Vector of molecules representing disconnected fragments
         */
        std::vector<Molecule> get_fragments() const;

        // Quantum chemistry compatibility
        /**
         * @brief Get atomic numbers as vector
         * @return Vector of atomic numbers in atom order
         */
        std::vector<int> atomic_numbers() const;
        
        /**
         * @brief Get element symbols as vector  
         * @return Vector of element symbols in atom order
         */
        std::vector<std::string> element_symbols() const;
        
        /**
         * @brief Get flattened geometry array for QC software
         * @return Flattened [x1,y1,z1,x2,y2,z2,...] coordinate array
         */
        xt::xarray<double> geometry_array() const;
        
        /**
         * @brief Get geometry as nested vector
         * @return Vector of [x,y,z] coordinate arrays
         */
        std::vector<std::array<double, 3>> geometry_vector() const;
        
        /**
         * @brief Get connectivity matrix for bond information
         * @return Matrix where [i][j] = bond order between atoms i and j
         */
        xt::xarray<double> connectivity_matrix() const;

        // File I/O and string representations
        /**
         * @brief Generate XYZ format string
         * @param include_comment Optional comment line
         * @return Multi-line XYZ format string
         */
        std::string to_xyz_string(const std::string& comment_line = "") const;
        
        /**
         * @brief Generate detailed molecular information string
         * @return Multi-line string with atoms, bonds, and properties
         */
        std::string to_detailed_string() const;
        
        /**
         * @brief Brief string representation
         * @return Single line molecular formula and basic info
         */
        std::string to_string() const;

        // Static construction methods
        /**
         * @brief Create molecule from XYZ format string
         * @param xyz_content Multi-line XYZ format string
         * @return Molecule object (bonds must be added separately)
         */
        static Molecule from_xyz_string(const std::string& xyz_content);
        
        /**
         * @brief Create simple linear molecule (for testing)
         * @param element_symbols Vector of element symbols
         * @param bond_length Distance between adjacent atoms
         * @return Linear molecule with single bonds
         */
        static Molecule create_linear_molecule(const std::vector<std::string>& element_symbols, 
                                             double bond_length = 1.5);

        // Molecular property storage
        /**
         * @brief Get stored molecular properties
         * @return Map of property names to values
         */
        const std::unordered_map<std::string, double>& properties() const { return molecular_properties_; }
        
        /**
         * @brief Get specific molecular property
         * @param property_name Name of property to retrieve
         * @return Property value if found, nullopt otherwise
         */
        std::optional<double> get_property(const std::string& property_name) const;
        
        /**
         * @brief Set molecular property
         * @param property_name Name of property
         * @param property_value Numerical value
         */
        void set_property(const std::string& property_name, double property_value);

        // Operators
        bool operator==(const Molecule& other) const;
        bool operator!=(const Molecule& other) const { return !(*this == other); }

    private:
        // Core data storage
        AtomContainer atoms_;                           ///< Unique atom instances
        BondContainer bonds_;                           ///< Bond instances with atom references
        int molecular_multiplicity_ = 1;               ///< Spin multiplicity (2S + 1)
        std::unordered_map<std::string, double> molecular_properties_;  ///< Quantum properties
        
        // Helper methods
        void validate_atom_index(AtomIndex index) const;
        void validate_bond_index(BondIndex index) const;
        bool bond_exists(AtomIndex first_index, AtomIndex second_index) const;
        void remove_bonds_involving_atom(AtomIndex atom_index);
        
        // Connectivity analysis helpers
        void depth_first_search_connectivity(AtomIndex current_atom, 
                                           std::vector<bool>& visited_atoms,
                                           std::vector<AtomIndex>& current_fragment) const;
    };

    // Implementation of simple inline methods

    inline const Atom& Molecule::atom(AtomIndex atom_index) const {
        validate_atom_index(atom_index);
        return *atoms_[atom_index];
    }

    inline Atom& Molecule::atom(AtomIndex atom_index) {
        validate_atom_index(atom_index);
        return *atoms_[atom_index];
    }

    inline const Bond& Molecule::bond(BondIndex bond_index) const {
        validate_bond_index(bond_index);
        return *bonds_[bond_index];
    }

    inline Bond& Molecule::bond(BondIndex bond_index) {
        validate_bond_index(bond_index);
        return *bonds_[bond_index];
    }

    inline AtomIndex Molecule::add_atom(int atomic_number, double x, double y, double z) {
        return add_atom(Atom(atomic_number, x, y, z));
    }

    inline std::optional<double> Molecule::get_property(const std::string& property_name) const {
        auto iterator = molecular_properties_.find(property_name);
        return (iterator != molecular_properties_.end()) ? 
               std::optional<double>(iterator->second) : std::nullopt;
    }

    inline void Molecule::set_property(const std::string& property_name, double property_value) {
        molecular_properties_[property_name] = property_value;
    }

    inline void Molecule::validate_atom_index(AtomIndex index) const {
        if (index >= atoms_.size()) {
            throw std::out_of_range("Atom index " + std::to_string(index) + 
                                  " out of range (molecule has " + std::to_string(atoms_.size()) + " atoms)");
        }
    }

    inline void Molecule::validate_bond_index(BondIndex index) const {
        if (index >= bonds_.size()) {
            throw std::out_of_range("Bond index " + std::to_string(index) + 
                                  " out of range (molecule has " + std::to_string(bonds_.size()) + " bonds)");
        }
    }

} // namespace libfrag

#endif // LIBFRAG_MOLECULE_HPP
