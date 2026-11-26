#!/usr/bin/env python3
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

"""Draws the AMCSET logo and writes every icon file the packages need.

The logo: an ion (amber) enters a crystal from the left, channels between two
rows of atoms, is deflected in two collisions and stops; the atoms it strikes
recoil (cyan). One geometry, in a 256 x 256 design square, is written as SVG
(the master, Linux's scalable icon) and drawn with Pillow for the bitmaps;
sizes of 32 pixels and less use a simplified drawing (no lattice, thicker
lines) so that they stay legible.

    python scripts/make_icons.py            # writes resources/icons/

Output (resources/icons/): amcset.svg; amcset-<n>.png for n = 16 ... 1024;
amcset.ico (16-256, for the Windows executables and installers); msix/ (the
tile and logo images of the MSIX package). Needs Pillow; the files are
committed, so a build does not run this script.
"""

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "resources" / "icons"

S = 256.0  # design square
RADIUS = 56.0  # corner radius of the tile
BG_TOP = (10, 27, 51)
BG_BOTTOM = (24, 74, 118)
MATERIAL = (255, 255, 255, 22)  # the crystal: a lighter region right of the surface
SURFACE_X = 62.0
ATOM = (128, 176, 222)
TRACK = (255, 184, 48)
RECOIL = (98, 214, 255)
WHITE = (255, 255, 255)

# Rows of atoms (y) and their columns (x): the ion channels between the rows
# at y = 92 and 128.
ROWS = [56.0, 92.0, 128.0, 164.0, 200.0]
COLS = [86.0, 122.0, 158.0, 194.0, 230.0]
# The ion: enters at the left, deflected at the two collision points, stops at the end.
ION = [(14.0, 110.0), (150.0, 110.0), (188.0, 162.0), (226.0, 178.0)]
# Recoils from the collision points.
RECOILS = [[(150.0, 110.0), (194.0, 70.0), (224.0, 58.0)], [(188.0, 162.0), (168.0, 214.0)]]
SMALL_RECOILS = [[(150.0, 110.0), (204.0, 62.0)], [(188.0, 162.0), (166.0, 220.0)]]


def lerp(a, b, t):
    return tuple(round(x + (y - x) * t) for x, y in zip(a, b))


def hexcolor(c):
    return "#%02x%02x%02x" % c[:3]


# --- SVG ----------------------------------------------------------------------


def svg():
    def path(points):
        return "M " + " L ".join(f"{x:g} {y:g}" for x, y in points)

    atoms = "\n".join(f'    <circle cx="{x:g}" cy="{y:g}" r="7"/>' for y in ROWS for x in COLS)
    recoils = "\n".join(f'    <path d="{path(r)}"/>' for r in RECOILS)
    m = MATERIAL
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 256 256" width="256" height="256">
  <title>AMCSET</title>
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0" stop-color="{hexcolor(BG_TOP)}"/>
      <stop offset="1" stop-color="{hexcolor(BG_BOTTOM)}"/>
    </linearGradient>
    <clipPath id="tile"><rect width="256" height="256" rx="{RADIUS:g}"/></clipPath>
    <filter id="glow" x="-20%" y="-20%" width="140%" height="140%">
      <feGaussianBlur stdDeviation="6"/>
    </filter>
  </defs>
  <rect width="256" height="256" rx="{RADIUS:g}" fill="url(#bg)"/>
  <g clip-path="url(#tile)">
    <rect x="{SURFACE_X:g}" y="0" width="{S - SURFACE_X:g}" height="256" fill="#ffffff" fill-opacity="{m[3] / 255:.3f}"/>
    <rect x="{SURFACE_X:g}" y="0" width="3" height="256" fill="#ffffff" fill-opacity="0.35"/>
  </g>
  <g fill="{hexcolor(ATOM)}" fill-opacity="0.55">
{atoms}
  </g>
  <g fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path d="{path(ION)}" stroke="{hexcolor(TRACK)}" stroke-width="26" stroke-opacity="0.45" filter="url(#glow)"/>
  </g>
  <g fill="none" stroke="{hexcolor(RECOIL)}" stroke-width="8" stroke-linecap="round" stroke-linejoin="round">
{recoils}
  </g>
  <path d="{path(ION)}" fill="none" stroke="{hexcolor(TRACK)}" stroke-width="13" stroke-linecap="round" stroke-linejoin="round"/>
  <g fill="#ffffff">
    <circle cx="{ION[1][0]:g}" cy="{ION[1][1]:g}" r="9"/>
    <circle cx="{ION[2][0]:g}" cy="{ION[2][1]:g}" r="9"/>
  </g>
  <circle cx="{ION[3][0]:g}" cy="{ION[3][1]:g}" r="12" fill="{hexcolor(TRACK)}" stroke="#ffffff" stroke-width="3"/>
</svg>
"""


# --- Bitmaps --------------------------------------------------------------------


def draw(size):
    """The logo at size x size pixels (drawn 4 times larger, then reduced)."""
    small = size <= 32
    k = 4
    n = size * k
    f = n / S

    def p(pt):
        return (pt[0] * f, pt[1] * f)

    def polyline(d, points, color, width):
        pts = [p(q) for q in points]
        d.line(pts, fill=color, width=max(1, round(width * f)), joint="curve")
        r = width * f / 2
        for x, y in (pts[0], pts[-1]):
            d.ellipse([x - r, y - r, x + r, y + r], fill=color)

    def dot(d, pt, r, color):
        x, y = p(pt)
        d.ellipse([x - r * f, y - r * f, x + r * f, y + r * f], fill=color)

    # Background gradient, clipped to the rounded tile.
    bg = Image.new("RGBA", (n, n))
    g = ImageDraw.Draw(bg)
    for y in range(n):
        g.line([(0, y), (n, y)], fill=lerp(BG_TOP, BG_BOTTOM, y / (n - 1)) + (255,))
    mask = Image.new("L", (n, n), 0)
    ImageDraw.Draw(mask).rounded_rectangle([0, 0, n - 1, n - 1], radius=RADIUS * f, fill=255)
    img = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    img.paste(bg, (0, 0), mask)

    # The crystal right of the surface.
    layer = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    d.rectangle([SURFACE_X * f, 0, n, n], fill=MATERIAL if not small else (255, 255, 255, 30))
    d.rectangle([SURFACE_X * f, 0, (SURFACE_X + (3 if not small else 6)) * f, n], fill=(255, 255, 255, 90))
    if not small:
        for y in ROWS:
            for x in COLS:
                dot(d, (x, y), 7, ATOM + (140,))
    clipped = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    clipped.paste(layer, (0, 0), mask)
    img = Image.alpha_composite(img, clipped)

    # Glow under the ion's track.
    if not small:
        glow = Image.new("RGBA", (n, n), (0, 0, 0, 0))
        polyline(ImageDraw.Draw(glow), ION, TRACK + (115,), 26)
        glow = glow.filter(ImageFilter.GaussianBlur(6 * f))
        img = Image.alpha_composite(img, glow)

    d = ImageDraw.Draw(img)
    for r in SMALL_RECOILS if small else RECOILS:
        polyline(d, r, RECOIL + (255,), 16 if small else 8)
    polyline(d, ION, TRACK + (255,), 24 if small else 13)
    if not small:
        dot(d, ION[1], 9, WHITE + (255,))
        dot(d, ION[2], 9, WHITE + (255,))
        dot(d, ION[3], 13.5, WHITE + (255,))
        dot(d, ION[3], 10.5, TRACK + (255,))
    else:
        dot(d, ION[3], 18, TRACK + (255,))
    return img.resize((size, size), Image.LANCZOS)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "amcset.svg").write_text(svg(), encoding="utf-8", newline="\n")
    images = {}
    for size in (16, 20, 22, 24, 30, 32, 36, 40, 44, 48, 50, 60, 64, 72, 80, 88, 96, 100, 128, 150, 256, 300,
                 512, 1024):
        images[size] = draw(size)
    for size in (16, 22, 24, 32, 48, 64, 128, 256, 512, 1024):
        images[size].save(OUT / f"amcset-{size}.png", optimize=True)
    ico = [16, 20, 24, 32, 40, 48, 64, 128, 256]
    images[256].save(OUT / "amcset.ico", sizes=[(s, s) for s in ico],
                     append_images=[images[s] for s in ico if s != 256])
    # MSIX: the app's logos at scale 100/200 and the taskbar sizes (targetsize).
    msix = OUT / "msix"
    msix.mkdir(exist_ok=True)
    for name, base in (("Square44x44Logo", 44), ("Square150x150Logo", 150), ("StoreLogo", 50)):
        for scale in (100, 200):
            images[base * scale // 100].save(msix / f"{name}.scale-{scale}.png", optimize=True)
    for t in (16, 24, 32, 48, 256):
        images[t].save(msix / f"Square44x44Logo.targetsize-{t}.png", optimize=True)
        images[t].save(msix / f"Square44x44Logo.targetsize-{t}_altform-unplated.png", optimize=True)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
