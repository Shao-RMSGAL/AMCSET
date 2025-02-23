import re
from pathlib import Path


def parse_value(value):
    # Remove uncertainty in parentheses and handle special cases
    if not value or value.isspace():
        return "0.0"  # Default value for empty fields

    # Extract the first number from range if present (e.g., [1.00784,1.00811])
    if value.startswith("["):
        match = re.search(r"\[([\d.]+)", value)
        if match:
            return match.group(1)
        return "0.0"

    # Remove uncertainty and '#' symbol
    match = re.search(r"([\d.]+)", value)
    if match:
        return match.group(1)
    return "0.0"


def generate_cpp_code(input_file):
    isotopes = []
    current_isotope = {}

    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_isotope:
                    isotopes.append(current_isotope)
                    current_isotope = {}
                continue

            if "=" in line:
                key, value = [x.strip() for x in line.split("=", 1)]
                current_isotope[key] = value

        if current_isotope:  # Don't forget the last isotope
            isotopes.append(current_isotope)

    # Generate C++ code
    cpp_code = """// Auto-generated code - DO NOT EDIT
#include <array>
#include <string>
#include "isotope_data.h"

namespace amcset {
namespace common {
    constexpr const IsotopeData::Isotope* IsotopeData::find_isotope(std::string_view atomic_symbol , int atomic_number) {
    constexpr std::array database = {
"""

    for isotope in isotopes:
        if not all(
            k in isotope for k in ["Atomic Symbol", "Mass Number", "Atomic Number"]
        ):
            continue

        symbol = isotope["Atomic Symbol"]
        mass_number = int(isotope["Mass Number"])

        cpp_code += f"""
        Isotope{{
        .atomic_number = {isotope['Atomic Number']},
        .atomic_symbol = "{symbol}",
        .mass_number = {mass_number},
        .relative_atomic_mass = {parse_value(isotope.get('Relative Atomic Mass', '0'))},
        .isotopic_composition = {parse_value(isotope.get('Isotopic Composition', '0'))},
        .standard_atomic_weight = {parse_value(isotope.get('Standard Atomic Weight', '0'))}
        }},
        """

    cpp_code += """
    }; 
    for (const auto& entry : database) {
        if (entry.atomic_symbol == atomic_symbol && entry.atomic_number == atomic_number) {
            return &entry;
        }
    }
    return nullptr;
}

} // namespace common
} // namespace amcset
"""
    return cpp_code


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python generate_isotope_data.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    cpp_code = generate_cpp_code(input_file)
    with open(output_file, "w") as f:
        f.write(cpp_code)
