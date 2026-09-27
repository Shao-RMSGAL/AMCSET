#!/usr/bin/env python3

from pathlib import Path
import sys

HEADER = """Copyright 2025, Texas A&M University

This file is part of AMCSET.

AMCSET is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
AMCSET. If not, see <https://www.gnu.org/licenses/>."""

EXTENSIONS = {".cpp", ".h",".py"}

EXCLUDED_DIRS = {
    ".git",
    "build",
    "generated",
}


def strip_comment_prefix(line: str) -> str:
    line = line.rstrip("\n")
    if line.startswith("//"):
        line = line[2:]
    elif line.startswith("#"):
        line = line[1:]
    if line.startswith(" "):
        line = line[1:]
    return line


def contains_header(text: str) -> bool:
    lines = text.split("\n")
    stripped_lines = [strip_comment_prefix(line) for line in lines]
    stripped_text = "\n".join(stripped_lines)
    return HEADER in stripped_text


failed = False

for path in Path(".").rglob("*"):
    if not path.is_file():
        continue

    if path.suffix.lower() not in EXTENSIONS:
        continue

    if any(part in EXCLUDED_DIRS for part in path.parts):
        continue

    text = path.read_text(encoding="utf-8", errors="replace")

    if not contains_header(text):
        print(f"ERROR: Missing header: {path}")
        failed = True

if failed:
    sys.exit(1)

print("All source files have the required header.")