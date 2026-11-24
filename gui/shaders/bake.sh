#!/bin/sh
# Re-compile the GUI shaders into shaders/baked/*.qsb (needs Qt Shader Tools).
#   QSB=/path/to/qsb gui/shaders/bake.sh
set -e
cd "$(dirname "$0")"
QSB=${QSB:-qsb}
for s in tracks.vert tracks.frag; do
  "$QSB" --glsl "100 es,120,150" --hlsl 50 --msl 12 -o "baked/$s.qsb" "$s"
done
