#pragma once
#ifndef LIBFRAG_ATOM_HPP
#define LIBFRAG_ATOM_HPP

#include <string>
#include <array>
#include <unordered_map>
#include <vector>
#include <optional>
#include <cmath>
#include <stdexcept>
#include "xtensor/xarray.hpp"

namespace libfrag {

    /**
     * @brief Comprehensive Atom class for quantum chemistry applications
     * 
     * Designed to interface with MolSSI Driver Interface, QCEngine, and QCElemental.
     * Stores atomic properties including coordinates, charge, mass, and quantum properties.
     */
    class Atom {
    public:
        // Constructors
        Atom() = default;
        
        /**
         * @brief Construct atom from atomic number and coordinates
         * @param atomic_number Atomic number (Z)
         * @param x X coordinate in Angstroms
         * @param y Y coordinate in Angstroms  
         * @param z Z coordinate in Angstroms
         */
        Atom(int atomic_number, double x, double y, double z);
        
        /**
         * @brief Construct atom from element symbol and coordinates
         * @param symbol Element symbol (e.g., "C", "N", "O")
         * @param x X coordinate in Angstroms
         * @param y Y coordinate in Angstroms
         * @param z Z coordinate in Angstroms
         */
        Atom(const std::string& symbol, double x, double y, double z);
        
        /**
         * @brief Construct atom from coordinates array
         * @param atomic_number Atomic number (Z)
         * @param coords Array of [x, y, z] coordinates in Angstroms
         */
        Atom(int atomic_number, const std::array<double, 3>& coords);
        
        /**
         * @brief Construct atom from xtensor array (for Python compatibility)
         * @param atomic_number Atomic number (Z)
         * @param coords xtensor array of coordinates
         */
        Atom(int atomic_number, const xt::xarray<double>& coords);

        // Accessors
        int atomic_number() const { return atomic_number_; }
        std::string symbol() const { return get_element_symbol(atomic_number_); }
        const std::array<double, 3>& coordinates() const { return coordinates_; }
        double x() const { return coordinates_[0]; }
        double y() const { return coordinates_[1]; }
        double z() const { return coordinates_[2]; }
        
        // Mass properties
        double atomic_mass() const { return get_atomic_mass(atomic_number_); }
        int mass_number() const { return mass_number_.value_or(get_most_common_isotope(atomic_number_)); }
        
        // Charge and electronic properties
        double formal_charge() const { return formal_charge_; }
        double partial_charge() const { return partial_charge_; }
        int multiplicity() const { return multiplicity_; }
        
        // Quantum properties
        const std::unordered_map<std::string, double>& properties() const { return properties_; }
        std::optional<double> get_property(const std::string& key) const;
        
        // Mutators
        void set_coordinates(double x, double y, double z);
        void set_coordinates(const std::array<double, 3>& coords);
        void set_coordinates(const xt::xarray<double>& coords);
        void set_formal_charge(double charge) { formal_charge_ = charge; }
        void set_partial_charge(double charge) { partial_charge_ = charge; }
        void set_multiplicity(int mult) { multiplicity_ = mult; }
        void set_mass_number(int mass_num) { mass_number_ = mass_num; }
        void set_property(const std::string& key, double value) { properties_[key] = value; }
        
        // Utility methods
        double distance_to(const Atom& other) const;
        std::array<double, 3> vector_to(const Atom& other) const;
        
        // Conversion methods for QC software compatibility
        xt::xarray<double> coordinates_array() const;
        std::vector<double> coordinates_vector() const;
        
        // String representations
        std::string to_string() const;
        std::string to_xyz_string() const;
        
        // Static utility methods
        static int symbol_to_atomic_number(const std::string& symbol);
        static std::string atomic_number_to_symbol(int atomic_number);
        static bool is_valid_atomic_number(int atomic_number);
        static bool is_valid_element_symbol(const std::string& symbol);
        
        // Operators
        bool operator==(const Atom& other) const;
        bool operator!=(const Atom& other) const { return !(*this == other); }

    private:
        int atomic_number_ = 1;  // Default to hydrogen
        std::array<double, 3> coordinates_ = {0.0, 0.0, 0.0};
        double formal_charge_ = 0.0;
        double partial_charge_ = 0.0;
        int multiplicity_ = 1;
        std::optional<int> mass_number_;
        std::unordered_map<std::string, double> properties_;
        
        // Static data for element properties
        static const std::unordered_map<std::string, int> symbol_to_z_;
        static const std::unordered_map<int, std::string> z_to_symbol_;
        static const std::unordered_map<int, double> atomic_masses_;
        static const std::unordered_map<int, int> most_common_isotopes_;
        
        // Helper methods
        static std::string get_element_symbol(int atomic_number);
        static double get_atomic_mass(int atomic_number);
        static int get_most_common_isotope(int atomic_number);
        static void initialize_element_data();
    };

    // Implementation of inline methods

    inline Atom::Atom(int atomic_number, double x, double y, double z)
        : atomic_number_(atomic_number), coordinates_({x, y, z}) {
        if (!is_valid_atomic_number(atomic_number)) {
            throw std::invalid_argument("Invalid atomic number: " + std::to_string(atomic_number));
        }
    }

    inline Atom::Atom(const std::string& symbol, double x, double y, double z)
        : coordinates_({x, y, z}) {
        atomic_number_ = symbol_to_atomic_number(symbol);
        if (atomic_number_ == 0) {
            throw std::invalid_argument("Invalid element symbol: " + symbol);
        }
    }

    inline Atom::Atom(int atomic_number, const std::array<double, 3>& coords)
        : atomic_number_(atomic_number), coordinates_(coords) {
        if (!is_valid_atomic_number(atomic_number)) {
            throw std::invalid_argument("Invalid atomic number: " + std::to_string(atomic_number));
        }
    }

    inline Atom::Atom(int atomic_number, const xt::xarray<double>& coords)
        : atomic_number_(atomic_number) {
        if (!is_valid_atomic_number(atomic_number)) {
            throw std::invalid_argument("Invalid atomic number: " + std::to_string(atomic_number));
        }
        if (coords.size() != 3) {
            throw std::invalid_argument("Coordinates array must have exactly 3 elements");
        }
        coordinates_[0] = coords(0);
        coordinates_[1] = coords(1);
        coordinates_[2] = coords(2);
    }

    inline std::optional<double> Atom::get_property(const std::string& key) const {
        auto it = properties_.find(key);
        return (it != properties_.end()) ? std::optional<double>(it->second) : std::nullopt;
    }

    inline void Atom::set_coordinates(double x, double y, double z) {
        coordinates_[0] = x;
        coordinates_[1] = y;
        coordinates_[2] = z;
    }

    inline void Atom::set_coordinates(const std::array<double, 3>& coords) {
        coordinates_ = coords;
    }

    inline void Atom::set_coordinates(const xt::xarray<double>& coords) {
        if (coords.size() != 3) {
            throw std::invalid_argument("Coordinates array must have exactly 3 elements");
        }
        coordinates_[0] = coords(0);
        coordinates_[1] = coords(1);
        coordinates_[2] = coords(2);
    }

    inline double Atom::distance_to(const Atom& other) const {
        double dx = coordinates_[0] - other.coordinates_[0];
        double dy = coordinates_[1] - other.coordinates_[1];
        double dz = coordinates_[2] - other.coordinates_[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    inline std::array<double, 3> Atom::vector_to(const Atom& other) const {
        return {
            other.coordinates_[0] - coordinates_[0],
            other.coordinates_[1] - coordinates_[1],
            other.coordinates_[2] - coordinates_[2]
        };
    }

    inline xt::xarray<double> Atom::coordinates_array() const {
        return xt::adapt(coordinates_);
    }

    inline std::vector<double> Atom::coordinates_vector() const {
        return {coordinates_[0], coordinates_[1], coordinates_[2]};
    }

    inline std::string Atom::to_string() const {
        return symbol() + " (" + std::to_string(coordinates_[0]) + ", " + 
               std::to_string(coordinates_[1]) + ", " + std::to_string(coordinates_[2]) + ")";
    }

    inline std::string Atom::to_xyz_string() const {
        char buffer[256];
        std::snprintf(buffer, sizeof(buffer), "%s %12.8f %12.8f %12.8f", 
                     symbol().c_str(), coordinates_[0], coordinates_[1], coordinates_[2]);
        return std::string(buffer);
    }

    inline bool Atom::operator==(const Atom& other) const {
        const double tolerance = 1e-10;
        return atomic_number_ == other.atomic_number_ &&
               std::abs(coordinates_[0] - other.coordinates_[0]) < tolerance &&
               std::abs(coordinates_[1] - other.coordinates_[1]) < tolerance &&
               std::abs(coordinates_[2] - other.coordinates_[2]) < tolerance &&
               std::abs(formal_charge_ - other.formal_charge_) < tolerance;
    }

    // Static data initialization (will be defined in implementation)
    inline const std::unordered_map<std::string, int> Atom::symbol_to_z_ = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8},
        {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},
        {"Cl", 17}, {"Ar", 18}, {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32},
        {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40},
        {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
        {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}
    };

    inline const std::unordered_map<int, std::string> Atom::z_to_symbol_ = {
        {1, "H"}, {2, "He"}, {3, "Li"}, {4, "Be"}, {5, "B"}, {6, "C"}, {7, "N"}, {8, "O"},
        {9, "F"}, {10, "Ne"}, {11, "Na"}, {12, "Mg"}, {13, "Al"}, {14, "Si"}, {15, "P"}, {16, "S"},
        {17, "Cl"}, {18, "Ar"}, {19, "K"}, {20, "Ca"}, {21, "Sc"}, {22, "Ti"}, {23, "V"}, {24, "Cr"},
        {25, "Mn"}, {26, "Fe"}, {27, "Co"}, {28, "Ni"}, {29, "Cu"}, {30, "Zn"}, {31, "Ga"}, {32, "Ge"},
        {33, "As"}, {34, "Se"}, {35, "Br"}, {36, "Kr"}, {37, "Rb"}, {38, "Sr"}, {39, "Y"}, {40, "Zr"},
        {41, "Nb"}, {42, "Mo"}, {43, "Tc"}, {44, "Ru"}, {45, "Rh"}, {46, "Pd"}, {47, "Ag"}, {48, "Cd"},
        {49, "In"}, {50, "Sn"}, {51, "Sb"}, {52, "Te"}, {53, "I"}, {54, "Xe"}
    };

    // Atomic masses (in atomic mass units)
    inline const std::unordered_map<int, double> Atom::atomic_masses_ = {
        {1, 1.008}, {2, 4.003}, {3, 6.941}, {4, 9.012}, {5, 10.811}, {6, 12.011}, {7, 14.007}, {8, 15.999},
        {9, 18.998}, {10, 20.180}, {11, 22.990}, {12, 24.305}, {13, 26.982}, {14, 28.086}, {15, 30.974}, {16, 32.066},
        {17, 35.453}, {18, 39.948}, {19, 39.098}, {20, 40.078}, {21, 44.956}, {22, 47.867}, {23, 50.942}, {24, 51.996},
        {25, 54.938}, {26, 55.845}, {27, 58.933}, {28, 58.693}, {29, 63.546}, {30, 65.380}
    };

    // Most common isotope mass numbers
    inline const std::unordered_map<int, int> Atom::most_common_isotopes_ = {
        {1, 1}, {2, 4}, {3, 7}, {4, 9}, {5, 11}, {6, 12}, {7, 14}, {8, 16},
        {9, 19}, {10, 20}, {11, 23}, {12, 24}, {13, 27}, {14, 28}, {15, 31}, {16, 32},
        {17, 35}, {18, 40}, {19, 39}, {20, 40}, {21, 45}, {22, 48}, {23, 51}, {24, 52},
        {25, 55}, {26, 56}, {27, 59}, {28, 58}, {29, 63}, {30, 64}
    };

    // Static method implementations
    inline int Atom::symbol_to_atomic_number(const std::string& symbol) {
        auto it = symbol_to_z_.find(symbol);
        return (it != symbol_to_z_.end()) ? it->second : 0;
    }

    inline std::string Atom::atomic_number_to_symbol(int atomic_number) {
        auto it = z_to_symbol_.find(atomic_number);
        return (it != z_to_symbol_.end()) ? it->second : "X";
    }

    inline bool Atom::is_valid_atomic_number(int atomic_number) {
        return atomic_number >= 1 && atomic_number <= 118;
    }

    inline bool Atom::is_valid_element_symbol(const std::string& symbol) {
        return symbol_to_z_.find(symbol) != symbol_to_z_.end();
    }

    inline std::string Atom::get_element_symbol(int atomic_number) {
        return atomic_number_to_symbol(atomic_number);
    }

    inline double Atom::get_atomic_mass(int atomic_number) {
        auto it = atomic_masses_.find(atomic_number);
        return (it != atomic_masses_.end()) ? it->second : 0.0;
    }

    inline int Atom::get_most_common_isotope(int atomic_number) {
        auto it = most_common_isotopes_.find(atomic_number);
        return (it != most_common_isotopes_.end()) ? it->second : atomic_number;
    }

} // namespace libfrag

#endif // LIBFRAG_ATOM_HPP
