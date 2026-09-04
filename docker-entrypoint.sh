#!/usr/bin/env bash
set -euo pipefail

# Run conan install now that we are at runtime and (hopefully) have a display.
if [[ -f "conanfile.py" ]]; then
    echo ">> Running conan install (build_type=Release)..."
    if ! conan install . --build=missing -s build_type=Release; then
        echo "!! conan install failed. Check the error above (display server? profile?)." >&2
        echo "!! Dropping you into a shell anyway so you can debug." >&2
    fi
else
    echo ">> No conanfile.py found in $(pwd); skipping conan install."
fi

exec "$@"
