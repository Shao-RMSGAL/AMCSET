#!/usr/bin/env bash
# Builds AMCSET-<version>-x86_64.AppImage (GUI + CLI + bundled SRIM tables).
#
# Needs a Linux machine (or WSL) with: GCC >= 13, CMake, Ninja, Qt 6 dev
# packages (Ubuntu: qt6-base-dev qt6-shadertools-dev libvulkan-dev libgl-dev
# libboost-dev libabsl-dev patchelf rsync libxcb-cursor0 file desktop-file-utils)
# and internet access for linuxdeploy (downloaded once into $WORK/tools).
# The binary is built with -march=x86-64-v2 (not native) so it runs on any
# recent x86-64 CPU. The AppImage carries the build machine's glibc: build on
# the oldest distribution you want to support.
#
#   tools/linux/make_appimage.sh [-o output_dir]
#
# Work folder (default ~/amcset-appimage; keep it off a synchronised folder):
#   AMCSET_APPIMAGE_WORK=/path tools/linux/make_appimage.sh
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
work="${AMCSET_APPIMAGE_WORK:-$HOME/amcset-appimage}"
out="$root/build"
while getopts "o:" opt; do case $opt in o) out="$OPTARG" ;; *) exit 2 ;; esac; done
version="$(sed -n 's/^project(amcset VERSION \([0-9.]*\)).*/\1/p' "$root/CMakeLists.txt")"
mkdir -p "$work/tools" "$out"

# 1. Source copy in a native folder (the repository may be on a slow or synchronised mount).
# (rsync keeps the timestamps of unchanged files, so a rerun rebuilds only what changed.)
mkdir -p "$work/src"
if git -C "$root" rev-parse --git-dir >/dev/null 2>&1; then
  git -C "$root" ls-files -z --cached --others --exclude-standard | rsync -a --from0 --files-from=- "$root/" "$work/src/"
else
  rsync -a --delete --exclude build --exclude .git "$root/" "$work/src/"
fi

# 2. Qt's private QRhi headers (package qt6-base-private-dev). Used from a
#    downloaded copy so that no root access is needed.
rhi="$(find /usr/include -path '*QtGui/*/QtGui/rhi/qrhi.h' 2>/dev/null | head -1 || true)"
if [ -z "$rhi" ]; then
  if [ ! -d "$work/qtprivate" ]; then
    (mkdir -p "$work/debs" && cd "$work/debs" && apt-get download qt6-base-private-dev \
      && dpkg -x qt6-base-private-dev_*.deb "$work/qtprivate")
  fi
  rhi="$(find "$work/qtprivate" -path '*QtGui/*/QtGui/rhi/qrhi.h' | head -1)"
fi
rhi_dir="${rhi%/rhi/qrhi.h}"

# 3. Build.
cmake -S "$work/src" -B "$work/build" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DAMCSET_MARCH=x86-64-v2 -DAMCSET_BUILD_TESTS=OFF -DAMCSET_BUILD_GUI=ON -DAMCSET_BUILD_CLI=ON \
  -DAMCSET_QT_RHI_INCLUDE_DIR="$rhi_dir"
ninja -C "$work/build"

# 4. AppDir.
appdir="$work/AppDir"
rm -rf "$appdir"
mkdir -p "$appdir/usr/bin" "$appdir/usr/share/amcset" \
         "$appdir/usr/share/applications" "$appdir/usr/share/icons/hicolor/scalable/apps" \
         "$appdir/usr/share/doc/amcset"
install -m 755 "$work/build/amcset-gui" "$work/build/amcset" "$appdir/usr/bin/"
cp -r "$work/src/data/srim2013" "$appdir/usr/share/amcset/"
cp "$work/src/README.md" "$work/src/LICENSE" "$work/src/CHANGELOG.md" "$appdir/usr/share/doc/amcset/" 2>/dev/null || true
cp -r "$work/src/docs" "$appdir/usr/share/doc/amcset/"
cp "$root/tools/linux/amcset.desktop" "$appdir/usr/share/applications/"
cp "$root/tools/linux/amcset.svg" "$appdir/usr/share/icons/hicolor/scalable/apps/amcset.svg"
# AppRun: the GUI by default; "AMCSET.AppImage amcset run scenario.json" etc. runs the CLI.
cat > "$work/AppRun" <<'EOF'
#!/bin/sh
HERE="$(dirname "$(readlink -f "$0")")"
export PATH="$HERE/usr/bin:$PATH"
if [ "${1:-}" = "amcset" ]; then shift; exec "$HERE/usr/bin/amcset" "$@"; fi
exec "$HERE/usr/bin/amcset-gui" "$@"
EOF
chmod +x "$work/AppRun"

# 5. linuxdeploy with the Qt plugin (downloaded once).
dl() { [ -x "$work/tools/$1" ] || { wget -q -O "$work/tools/$1" "$2" && chmod +x "$work/tools/$1"; }; }
dl linuxdeploy-x86_64.AppImage \
   https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
dl linuxdeploy-plugin-qt-x86_64.AppImage \
   https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/continuous/linuxdeploy-plugin-qt-x86_64.AppImage
export APPIMAGE_EXTRACT_AND_RUN=1 NO_STRIP=1 PATH="$work/tools:$PATH"
export QMAKE="${QMAKE:-$(command -v qmake6 || echo /usr/lib/qt6/bin/qmake6)}"
export LINUXDEPLOY_OUTPUT_VERSION="$version"
cd "$work"
rm -f ./*.AppImage
"$work/tools/linuxdeploy-x86_64.AppImage" --appdir "$appdir" \
  --executable "$appdir/usr/bin/amcset-gui" --executable "$appdir/usr/bin/amcset" \
  --desktop-file "$appdir/usr/share/applications/amcset.desktop" \
  --icon-file "$appdir/usr/share/icons/hicolor/scalable/apps/amcset.svg" \
  --custom-apprun "$work/AppRun" --plugin qt --output appimage

img="$(ls -1 "$work"/*.AppImage | head -1)"
cp "$img" "$out/AMCSET-$version-x86_64.AppImage"
echo "Built $out/AMCSET-$version-x86_64.AppImage ($(du -h "$out/AMCSET-$version-x86_64.AppImage" | cut -f1))"
