#!/usr/bin/env bash
#
# One-off script: add the AMCSET license header to every existing
# .h / .cpp / .py file under visualization/, any src/ or include/
# directory, or test/common/ that doesn't already have it.
#
# Run from the repo root:
#   chmod +x add-headers.sh
#   ./add-headers.sh
#
# Then review with `git diff` and commit as usual.

set -euo pipefail

read -r -d '' HEADER_CPP <<'EOF' || true
// Copyright 2025, Texas A&M University
//
// This file is part of AMCSET.
//
// AMCSET is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.
EOF

read -r -d '' HEADER_PY <<'EOF' || true
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
EOF

path_requires_header() {
    local f="$1"
    case "$f" in
        visualization/*) return 0 ;;
        test/common/*)   return 0 ;;
        */src/*|src/*)         return 0 ;;
        */include/*|include/*) return 0 ;;
        *) return 1 ;;
    esac
}

updated=0
skipped=0

# Walk all tracked files (respects .gitignore); fall back to `find`
# if you're not in a git repo or want untracked files included too.
while IFS= read -r -d '' file; do
    [ -f "$file" ] || continue

    case "$file" in
        *.h|*.cpp) header="$HEADER_CPP" ;;
        *.py)      header="$HEADER_PY" ;;
        *) continue ;;
    esac

    path_requires_header "$file" || continue

    header_lines=$(printf '%s\n' "$header" | wc -l)
    # Strip \r so CRLF-vs-LF differences don't cause false mismatches.
    actual=$(head -n "$header_lines" "$file" | tr -d '\r')
    normalized_header=$(printf '%s\n' "$header" | tr -d '\r')

    if [ "$actual" = "$normalized_header" ]; then
        skipped=$((skipped + 1))
        continue
    fi

    # Belt-and-suspenders: also skip if the header text appears anywhere
    # in the first ~20 lines (handles a stray blank line/BOM before it).
    if head -n 20 "$file" | tr -d '\r' | grep -qF "Copyright 2025, Texas A&M University"; then
        skipped=$((skipped + 1))
        continue
    fi

    tmp=$(mktemp)
    {
        printf '%s\n' "$header"
        if [ -s "$file" ] && [ "$(head -c1 "$file")" != $'\n' ]; then
            printf '\n'
        fi
        cat "$file"
    } > "$tmp"
    mv "$tmp" "$file"
    updated=$((updated + 1))
    echo "Added header: $file"
done < <(git ls-files -z -- '*.h' '*.cpp' '*.py')

echo ""
echo "Done. Headers added: $updated, already present: $skipped."
echo "Review the changes with 'git diff', then 'git add' and commit."