# Copyright 2025, Texas A&M University
#
# This file is part of AMCSET.
#
# AMCSET is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# AMCSET. If not, see <https://www.gnu.org/licenses/>.

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
#pragma once
#define GENERATE_ISOTOPE_DATA \\
"""

    for isotope in isotopes:
        if not all(
            k in isotope for k in ["Atomic Symbol", "Mass Number", "Atomic Number"]
        ):
            continue

        symbol = isotope["Atomic Symbol"]
        mass_number = int(isotope["Mass Number"])

        cpp_code += f"""Isotope{{ \\
        .atomic_number = {isotope['Atomic Number']}, \\
        .atomic_symbol = "{symbol}", \\
        .mass_number = {mass_number}, \\
        .relative_atomic_mass = {parse_value(isotope.get('Relative Atomic Mass', '0'))}, \\
        .isotopic_composition = {parse_value(isotope.get('Isotopic Composition', '0'))}, \\
        .standard_atomic_weight = {parse_value(isotope.get('Standard Atomic Weight', '0'))} \\
        }}, \\
        """

    cpp_code += """
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
