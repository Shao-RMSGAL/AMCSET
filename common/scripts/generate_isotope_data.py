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
#include <map>
#include <string>
#include "isotope.h"

namespace chemistry {

struct Isotope {
    int atomic_number;
    std::string atomic_symbol;
    int mass_number;
    double relative_atomic_mass;
    double isotopic_composition;
    double standard_atomic_weight;
};

std::map<std::string, std::map<int, Isotope>> createIsotopeDatabase() {
    std::map<std::string, std::map<int, Isotope>> database;
"""

    for isotope in isotopes:
        if not all(
            k in isotope for k in ["Atomic Symbol", "Mass Number", "Atomic Number"]
        ):
            continue

        symbol = isotope["Atomic Symbol"]
        mass_number = int(isotope["Mass Number"])

        cpp_code += f"""    database["{symbol}"][{mass_number}] = Isotope{{
        .atomic_number = {isotope['Atomic Number']},
        .atomic_symbol = "{symbol}",
        .mass_number = {mass_number},
        .relative_atomic_mass = {parse_value(isotope.get('Relative Atomic Mass', '0'))},
        .isotopic_composition = {parse_value(isotope.get('Isotopic Composition', '0'))},
        .standard_atomic_weight = {parse_value(isotope.get('Standard Atomic Weight', '0'))}
    }};\n"""

    cpp_code += """
    return database;
}

} // namespace chemistry
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

#  # Add the custom command to generate the isotope database
#  add_custom_command(
#      OUTPUT ${CMAKE_BINARY_DIR}/generated/isotope_database.cpp
#      COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_isotope_data.py
#              ${CMAKE_SOURCE_DIR}/assets/isotope_data.txt
#              ${CMAKE_BINARY_DIR}/generated/isotope_database.cpp
#      DEPENDS ${CMAKE_SOURCE_DIR}/assets/isotope_data.txt
#              ${CMAKE_SOURCE_DIR}/scripts/generate_isotope_data.py
#      COMMENT "Generating isotope database"
#  )

#  # Create a library that includes the generated code
#  add_library(isotope_database
#      ${CMAKE_BINARY_DIR}/generated/isotope_database.cpp
#      include/isotope.h
#  )

#  target_include_directories(isotope_database
#      PUBLIC
#          ${CMAKE_SOURCE_DIR}/include
#          ${CMAKE_BINARY_DIR}/generated
#  )
