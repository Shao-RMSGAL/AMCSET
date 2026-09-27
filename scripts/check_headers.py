#!/usr/bin/env python3

from pathlib import Path
import sys

HEADER = "Copyright (c) 2026 Nathaniel Thomas"

EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx"}

EXCLUDED_DIRS = {
    ".git",
    "build",
    "third_party",
    "external",
}

failed = False

for path in Path(".").rglob("*"):
    if not path.is_file():
        continue

    if path.suffix.lower() not in EXTENSIONS:
        continue

    if any(part in EXCLUDED_DIRS for part in path.parts):
        continue

    text = path.read_text(encoding="utf-8", errors="replace")

    if HEADER not in text:
        print(f"ERROR: Missing header: {path}")
        failed = True

if failed:
    sys.exit(1)

print("All source files have the required header.")