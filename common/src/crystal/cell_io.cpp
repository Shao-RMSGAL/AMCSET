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

#include "crystal/cell_io.h"

#include "element_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace crystal {

using namespace common;

namespace {

std::string lower(std::string_view s) {
  std::string o(s);
  for (auto &c : o) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return o;
}

std::string trim(std::string_view s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return std::string(s.substr(a, b - a));
}

std::vector<std::string> split_ws(std::string_view s) {
  std::vector<std::string> out;
  std::istringstream in{std::string(s)};
  std::string t;
  while (in >> t) out.push_back(t);
  return out;
}

std::vector<std::string> lines_of(std::string_view text) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : text) {
    if (c == '\n') {
      if (!cur.empty() && cur.back() == '\r') cur.pop_back();
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  if (!cur.empty()) {
    if (cur.back() == '\r') cur.pop_back();
    out.push_back(cur);
  }
  return out;
}

//! A number; \p what names it in the message. Accepts CIF's standard
//! uncertainty in parentheses ("5.4310(2)").
double to_number(std::string_view text, const std::string &what) {
  std::string t = trim(text);
  if (const auto p = t.find('('); p != std::string::npos) t = t.substr(0, p);
  if (t.empty()) throw std::invalid_argument(what + ": missing number");
  std::size_t used = 0;
  double v = 0.0;
  try {
    v = std::stod(t, &used);
  } catch (const std::exception &) {
    throw std::invalid_argument(what + ": '" + std::string(text) + "' is not a number");
  }
  if (used != t.size() || !std::isfinite(v)) {
    throw std::invalid_argument(what + ": '" + std::string(text) + "' is not a number");
  }
  return v;
}

std::string at_line(const std::string &source, std::size_t line) {
  return source + ", line " + std::to_string(line + 1);
}

//! Token \p i of \p t. \throws std::invalid_argument if the line is too short.
const std::string &token(const std::vector<std::string> &t, std::size_t i, const std::string &where) {
  if (i >= t.size()) throw std::invalid_argument(where + ": expected at least " + std::to_string(i + 1) + " value(s)");
  return t[i];
}

//! A position of one element read from a file, before merging into sites.
struct RawAtom {
  std::string label;
  unsigned z = 0;
  double fraction = 1.0;
  Vec3 position;               // fractional
  SiteKind kind = SiteKind::lattice;
  std::size_t row = 0;         // source row (CIF): images of one row are deduplicated
};

/*!
 * Merge atoms into sites: images of the same source row closer than tol are
 * duplicates (special positions); different rows at the same position are
 * one site with several occupants (fractions of the same element add up).
 */
std::vector<Site> merge_atoms(const Lattice &lattice, const std::vector<RawAtom> &atoms, double tol_A,
                              const std::string &source) {
  UnitCell probe;
  probe.lattice = lattice;
  std::vector<Site> sites;
  std::vector<std::vector<std::size_t>> rows_at; // rows already on each site
  for (const auto &a : atoms) {
    const Vec3 p = wrap_fractional(a.position);
    std::size_t found = sites.size();
    for (std::size_t i = 0; i < sites.size(); ++i) {
      if (probe.distance_A(sites[i].position, p) < tol_A) {
        found = i;
        break;
      }
    }
    if (found == sites.size()) {
      Site s;
      s.label = a.label;
      s.position = p;
      s.kind = a.kind;
      s.occupants.push_back({a.z, a.fraction});
      sites.push_back(std::move(s));
      rows_at.push_back({a.row});
      continue;
    }
    auto &rows = rows_at[found];
    if (std::find(rows.begin(), rows.end(), a.row) != rows.end()) continue; // an image of the same row
    rows.push_back(a.row);
    auto &occ = sites[found].occupants;
    auto it = std::find_if(occ.begin(), occ.end(), [&](const Occupant &o) { return o.z == a.z; });
    if (it == occ.end()) {
      occ.push_back({a.z, a.fraction});
    } else {
      it->fraction += a.fraction;
    }
    if (a.kind == SiteKind::interstitial) sites[found].kind = SiteKind::interstitial;
    if (sites[found].occupancy() > 1.0 + 1e-3) {
      throw std::invalid_argument(source + ": the occupancies at site " + sites[found].label + " add up to more than 1");
    }
  }
  // Rounded file values (0.333 + 0.667) may add up to slightly above 1.
  for (auto &s : sites) {
    const double occ = s.occupancy();
    if (occ > 1.0) {
      for (auto &o : s.occupants) o.fraction /= occ;
    }
  }
  return sites;
}

std::string symbol(unsigned z) { return std::string(Element::by_z(z).symbol()); }

//! Fixed notation with \p decimals digits (no "-0.000").
std::string fmt(double v, int decimals = 10) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*f", decimals, std::abs(v) < 5e-13 ? 0.0 : v);
  return buf;
}

//! The shortest of fixed and exponent notation with \p digits significant digits.
std::string fmt_general(double v, int digits) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*g", digits, v);
  return buf;
}

//! Elements in order of first appearance, and whether every site is fully
//! occupied by one element (what POSCAR, XYZ and LAMMPS data can store).
std::vector<unsigned> species_order(const UnitCell &cell, const char *format) {
  std::vector<unsigned> order;
  for (const auto &s : cell.sites) {
    if (s.occupants.empty()) continue;
    if (s.occupants.size() != 1 || std::abs(s.occupants.front().fraction - 1.0) > 1e-9) {
      throw std::invalid_argument(std::string(format) + " cannot store the partly occupied site " + s.label +
                                  "; write a CIF instead");
    }
    const unsigned z = s.occupants.front().z;
    if (std::find(order.begin(), order.end(), z) == order.end()) order.push_back(z);
  }
  if (order.empty()) throw std::invalid_argument(std::string(format) + ": the cell has no atoms");
  return order;
}

// ---------------------------------------------------------------------------
// CIF
// ---------------------------------------------------------------------------

struct CifToken {
  std::string text;
  bool quoted = false;
  std::size_t line = 0;
};

std::vector<CifToken> cif_tokens(std::string_view text, const std::string &source) {
  std::vector<CifToken> out;
  const auto lines = lines_of(text);
  for (std::size_t ln = 0; ln < lines.size(); ++ln) {
    const std::string &l = lines[ln];
    if (!l.empty() && l[0] == ';') {
      // Text field up to the next line that starts with ';'.
      std::string field = l.substr(1);
      std::size_t end = ln + 1;
      while (end < lines.size() && (lines[end].empty() || lines[end][0] != ';')) {
        field += "\n" + lines[end];
        ++end;
      }
      if (end >= lines.size()) throw std::invalid_argument(at_line(source, ln) + ": unterminated ; text field");
      out.push_back({trim(field), true, ln});
      ln = end; // the closing ';' line (anything after it on that line is ignored)
      continue;
    }
    std::size_t i = 0;
    while (i < l.size()) {
      const char c = l[i];
      if (std::isspace(static_cast<unsigned char>(c))) {
        ++i;
        continue;
      }
      if (c == '#') break;
      if (c == '\'' || c == '"') {
        // A quote closes only when followed by white space or the line end.
        std::size_t j = i + 1;
        while (j < l.size() && !(l[j] == c && (j + 1 == l.size() || std::isspace(static_cast<unsigned char>(l[j + 1])))))
          ++j;
        if (j >= l.size()) throw std::invalid_argument(at_line(source, ln) + ": unterminated quoted string");
        out.push_back({l.substr(i + 1, j - i - 1), true, ln});
        i = j + 1;
        continue;
      }
      std::size_t j = i;
      while (j < l.size() && !std::isspace(static_cast<unsigned char>(l[j]))) ++j;
      out.push_back({l.substr(i, j - i), false, ln});
      i = j;
    }
  }
  return out;
}

struct CifLoop {
  std::vector<std::string> tags; // lower case
  std::vector<std::vector<CifToken>> rows;
  int column(std::string_view tag) const {
    for (std::size_t i = 0; i < tags.size(); ++i) {
      if (tags[i] == tag) return static_cast<int>(i);
    }
    return -1;
  }
};

struct CifBlock {
  std::string name;
  std::map<std::string, CifToken> items; // lower-case tag -> value
  std::vector<CifLoop> loops;
  const CifLoop *loop_with(std::string_view tag) const {
    for (const auto &l : loops) {
      if (l.column(tag) >= 0) return &l;
    }
    return nullptr;
  }
};

bool is_reserved(const CifToken &t) {
  if (t.quoted) return false;
  const std::string l = lower(t.text);
  return l == "loop_" || l.rfind("data_", 0) == 0 || l == "global_" || l.rfind("save_", 0) == 0 || l == "stop_";
}

CifBlock parse_cif_block(std::string_view text, const std::string &source) {
  const auto tok = cif_tokens(text, source);
  CifBlock b;
  std::size_t i = 0;
  bool in_block = false;
  while (i < tok.size()) {
    const CifToken &t = tok[i];
    const std::string lt = lower(t.text);
    if (!t.quoted && lt.rfind("data_", 0) == 0) {
      if (in_block) break; // only the first block
      in_block = true;
      b.name = t.text.substr(5);
      ++i;
      continue;
    }
    if (!in_block) {
      ++i; // anything before the first data_ (should be nothing)
      continue;
    }
    if (!t.quoted && lt == "loop_") {
      CifLoop loop;
      ++i;
      while (i < tok.size() && !tok[i].quoted && !tok[i].text.empty() && tok[i].text[0] == '_') {
        loop.tags.push_back(lower(tok[i].text));
        ++i;
      }
      if (loop.tags.empty()) throw std::invalid_argument(at_line(source, t.line) + ": loop_ without tags");
      std::vector<CifToken> values;
      while (i < tok.size() && !is_reserved(tok[i]) && (tok[i].quoted || tok[i].text.empty() || tok[i].text[0] != '_')) {
        values.push_back(tok[i]);
        ++i;
      }
      if (values.size() % loop.tags.size() != 0) {
        throw std::invalid_argument(at_line(source, t.line) + ": a loop with " + std::to_string(loop.tags.size()) +
                                    " tags has " + std::to_string(values.size()) + " values (not a multiple)");
      }
      for (std::size_t k = 0; k < values.size(); k += loop.tags.size()) {
        loop.rows.emplace_back(values.begin() + static_cast<std::ptrdiff_t>(k),
                               values.begin() + static_cast<std::ptrdiff_t>(k + loop.tags.size()));
      }
      b.loops.push_back(std::move(loop));
      continue;
    }
    if (!t.quoted && !t.text.empty() && t.text[0] == '_') {
      if (i + 1 >= tok.size() || is_reserved(tok[i + 1]) ||
          (!tok[i + 1].quoted && !tok[i + 1].text.empty() && tok[i + 1].text[0] == '_')) {
        throw std::invalid_argument(at_line(source, t.line) + ": " + t.text + " has no value");
      }
      b.items[lt] = tok[i + 1];
      i += 2;
      continue;
    }
    ++i; // stray value: ignore
  }
  if (!in_block) throw std::invalid_argument(source + ": no data_ block (not a CIF?)");
  return b;
}

bool cif_missing(const CifToken &t) { return !t.quoted && (t.text == "." || t.text == "?"); }

double cif_number(const CifBlock &b, const std::string &tag, const std::string &source, std::optional<double> def) {
  const auto it = b.items.find(tag);
  if (it == b.items.end() || cif_missing(it->second)) {
    if (def) return *def;
    throw std::invalid_argument(source + ": " + tag + " is missing");
  }
  return to_number(it->second.text, at_line(source, it->second.line) + ": " + tag);
}

UnitCell parse_cif(std::string_view text, const CellReadOptions &opt, const std::string &source) {
  const CifBlock b = parse_cif_block(text, source);
  UnitCell cell;
  cell.name = b.name;
  if (const auto it = b.items.find("_chemical_formula_sum"); it != b.items.end() && !cif_missing(it->second)) {
    cell.name = trim(it->second.text) + (b.name.empty() ? "" : " (" + b.name + ")");
  }
  cell.lattice = Lattice::from_parameters(
      cif_number(b, "_cell_length_a", source, std::nullopt) * angstrom,
      cif_number(b, "_cell_length_b", source, std::nullopt) * angstrom,
      cif_number(b, "_cell_length_c", source, std::nullopt) * angstrom,
      cif_number(b, "_cell_angle_alpha", source, 90.0) * degree,
      cif_number(b, "_cell_angle_beta", source, 90.0) * degree,
      cif_number(b, "_cell_angle_gamma", source, 90.0) * degree);

  // Symmetry operations.
  std::vector<std::pair<Mat3, Vec3>> ops;
  for (const char *tag : {"_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"}) {
    if (const CifLoop *l = b.loop_with(tag)) {
      const int c = l->column(tag);
      for (const auto &row : l->rows) {
        try {
          ops.push_back(parse_symmetry_operation(row[static_cast<std::size_t>(c)].text));
        } catch (const std::invalid_argument &e) {
          throw std::invalid_argument(at_line(source, row[static_cast<std::size_t>(c)].line) + ": " + e.what());
        }
      }
      break;
    }
    if (const auto it = b.items.find(tag); it != b.items.end()) {
      ops.push_back(parse_symmetry_operation(it->second.text));
      break;
    }
  }
  if (ops.empty()) {
    std::string group;
    for (const char *tag : {"_symmetry_space_group_name_h-m", "_space_group_name_h-m_alt", "_space_group_name_hall",
                            "_symmetry_space_group_name_hall"}) {
      if (const auto it = b.items.find(tag); it != b.items.end() && !cif_missing(it->second)) group = it->second.text;
    }
    std::string g;
    for (char ch : group) {
      if (!std::isspace(static_cast<unsigned char>(ch))) g += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    for (const char *tag : {"_space_group_it_number", "_symmetry_int_tables_number"}) {
      const auto it = b.items.find(tag);
      if (it == b.items.end() || cif_missing(it->second)) continue;
      if (to_number(it->second.text, at_line(source, it->second.line) + ": " + tag) != 1.0 && g.empty()) {
        g = "number " + trim(it->second.text);
        group = g;
      }
    }
    if (!g.empty() && g != "p1") {
      throw std::invalid_argument(source + ": space group '" + group +
                                  "' without its symmetry operations (_symmetry_equiv_pos_as_xyz or "
                                  "_space_group_symop_operation_xyz); export the CIF with them, or in P 1");
    }
    ops.push_back({Mat3::identity(), Vec3{}});
  }

  // Atom sites.
  const CifLoop *sites = b.loop_with("_atom_site_fract_x");
  bool cartesian = false;
  if (!sites) {
    sites = b.loop_with("_atom_site_cartn_x");
    cartesian = sites != nullptr;
  }
  if (!sites) throw std::invalid_argument(source + ": no atom sites (_atom_site_fract_x or _atom_site_Cartn_x)");
  const std::string px = cartesian ? "_atom_site_cartn_" : "_atom_site_fract_";
  const int cx = sites->column(px + "x"), cy = sites->column(px + "y"), cz = sites->column(px + "z");
  if (cy < 0 || cz < 0) throw std::invalid_argument(source + ": atom sites without all three coordinates");
  const int c_type = sites->column("_atom_site_type_symbol");
  const int c_label = sites->column("_atom_site_label");
  const int c_occ = sites->column("_atom_site_occupancy");
  const int c_kind = sites->column("_amcset_site_kind");
  std::vector<RawAtom> atoms;
  for (std::size_t r = 0; r < sites->rows.size(); ++r) {
    const auto &row = sites->rows[r];
    const auto cell_at = [&](int c) -> const CifToken & { return row[static_cast<std::size_t>(c)]; };
    const std::string where = at_line(source, row.front().line);
    std::string label = c_label >= 0 ? cell_at(c_label).text : "";
    std::optional<unsigned> z;
    if (c_type >= 0 && !cif_missing(cell_at(c_type))) z = element_of_symbol(cell_at(c_type).text);
    if (!z && !label.empty()) z = element_of_symbol(label);
    if (!z) throw std::invalid_argument(where + ": no element for atom site '" + label + "'");
    if (label.empty()) label = symbol(*z) + std::to_string(r + 1);
    double occ = 1.0;
    if (c_occ >= 0 && !cif_missing(cell_at(c_occ))) occ = to_number(cell_at(c_occ).text, where + ": occupancy");
    if (!(occ > 0.0) || occ > 1.0 + 1e-6) {
      throw std::invalid_argument(where + ": occupancy of " + label + " must lie in (0, 1]");
    }
    occ = std::min(occ, 1.0);
    Vec3 p{to_number(cell_at(cx).text, where + ": x"), to_number(cell_at(cy).text, where + ": y"),
           to_number(cell_at(cz).text, where + ": z")};
    if (cartesian) p = cell.lattice.to_fractional(p);
    SiteKind kind = SiteKind::lattice;
    if (c_kind >= 0 && !cif_missing(cell_at(c_kind))) kind = parse_site_kind(cell_at(c_kind).text);
    for (const auto &[m, t] : ops) atoms.push_back({label, *z, occ, m * p + t, kind, r});
  }
  cell.sites = merge_atoms(cell.lattice, atoms, opt.merge_distance_A, source);
  return cell;
}

// ---------------------------------------------------------------------------
// POSCAR
// ---------------------------------------------------------------------------

bool starts_alpha(const std::string &s) { return !s.empty() && std::isalpha(static_cast<unsigned char>(s[0])); }

UnitCell parse_poscar(std::string_view text, const CellReadOptions &opt, const std::string &source) {
  const auto lines = lines_of(text);
  std::size_t ln = 0;
  const auto next = [&]() -> const std::string & {
    if (ln >= lines.size()) throw std::invalid_argument(source + ": ends early (line " + std::to_string(ln + 1) + ")");
    return lines[ln++];
  };
  UnitCell cell;
  cell.name = trim(next());
  const auto scale_line = split_ws(next());
  const double scale = to_number(token(scale_line, 0, at_line(source, ln - 1) + ": scale"),
                                 at_line(source, ln - 1) + ": scale");
  Mat3 m;
  for (int i = 0; i < 3; ++i) {
    const auto t = split_ws(next());
    if (t.size() < 3) throw std::invalid_argument(at_line(source, ln - 1) + ": expected a lattice vector");
    for (int j = 0; j < 3; ++j) m[i][j] = to_number(t[static_cast<std::size_t>(j)], at_line(source, ln - 1) + ": lattice");
  }
  double factor = scale;
  if (scale < 0.0) factor = std::cbrt(-scale / std::abs(m.determinant())); // VASP: negative = volume
  if (!(factor > 0.0)) throw std::invalid_argument(at_line(source, 1) + ": the scale must not be 0");
  for (int i = 0; i < 3; ++i) m[i] = m[i] * factor;
  cell.lattice = Lattice::from_vectors_A(m);
  auto t = split_ws(next());
  std::vector<unsigned> species;
  if (!t.empty() && starts_alpha(t[0])) {
    for (const auto &s : t) {
      // "Fe_pv", "Fe/abc..." (POTCAR names) -> Fe
      const std::string base = s.substr(0, s.find_first_of("_/"));
      const auto z = Element::find_z(base);
      if (!z) throw std::invalid_argument(at_line(source, ln - 1) + ": unknown element '" + s + "'");
      species.push_back(*z);
    }
    t = split_ws(next());
  } else {
    species = opt.elements;
    if (species.empty()) {
      throw std::invalid_argument(source + ": a VASP 4 POSCAR without a species line; give the elements in the "
                                           "order of the counts");
    }
  }
  std::vector<std::size_t> counts;
  for (const auto &c : t) {
    const double v = to_number(c, at_line(source, ln - 1) + ": counts");
    if (v < 0 || v != std::floor(v)) throw std::invalid_argument(at_line(source, ln - 1) + ": counts must be whole numbers");
    counts.push_back(static_cast<std::size_t>(v));
  }
  if (counts.size() != species.size()) {
    throw std::invalid_argument(at_line(source, ln - 1) + ": " + std::to_string(counts.size()) + " counts for " +
                                std::to_string(species.size()) + " elements");
  }
  std::string mode = trim(next());
  if (!mode.empty() && (mode[0] == 'S' || mode[0] == 's')) mode = trim(next()); // Selective dynamics
  const bool cart = !mode.empty() && (mode[0] == 'C' || mode[0] == 'c' || mode[0] == 'K' || mode[0] == 'k');
  if (!cart && !(mode.empty() || mode[0] == 'D' || mode[0] == 'd')) {
    throw std::invalid_argument(at_line(source, ln - 1) + ": expected Direct or Cartesian, got '" + mode + "'");
  }
  std::vector<RawAtom> atoms;
  std::size_t row = 0;
  for (std::size_t s = 0; s < species.size(); ++s) {
    for (std::size_t k = 0; k < counts[s]; ++k) {
      const auto c = split_ws(next());
      if (c.size() < 3) throw std::invalid_argument(at_line(source, ln - 1) + ": expected three coordinates");
      Vec3 p{to_number(c[0], at_line(source, ln - 1)), to_number(c[1], at_line(source, ln - 1)),
             to_number(c[2], at_line(source, ln - 1))};
      if (cart) p = cell.lattice.to_fractional(p * factor);
      atoms.push_back({symbol(species[s]) + std::to_string(k + 1), species[s], 1.0, p, SiteKind::lattice, row++});
    }
  }
  if (atoms.empty()) throw std::invalid_argument(source + ": no atoms");
  cell.sites = merge_atoms(cell.lattice, atoms, opt.merge_distance_A, source);
  return cell;
}

// ---------------------------------------------------------------------------
// LAMMPS data
// ---------------------------------------------------------------------------

bool numeric_start(const std::string &s) {
  const std::string t = trim(s);
  return !t.empty() && (std::isdigit(static_cast<unsigned char>(t[0])) || t[0] == '-' || t[0] == '+' || t[0] == '.');
}

std::string strip_comment(const std::string &s) { return trim(s.substr(0, s.find('#'))); }

UnitCell parse_lammps_data(std::string_view text, const CellReadOptions &opt, const std::string &source) {
  const auto lines = lines_of(text);
  UnitCell cell;
  cell.name = lines.empty() ? "" : trim(lines[0]);
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
  bool have[3] = {false, false, false};
  double xy = 0, xz = 0, yz = 0;
  std::size_t n_atoms = 0, n_types = 0;
  std::size_t ln = 1;
  // Header: lines with numbers first, up to the first section keyword.
  for (; ln < lines.size(); ++ln) {
    const std::string l = strip_comment(lines[ln]);
    if (l.empty()) continue;
    if (!numeric_start(l)) break;
    const auto t = split_ws(l);
    const std::string ll = lower(l);
    const auto ends = [&](std::string_view s) {
      return ll.size() >= s.size() && ll.compare(ll.size() - s.size(), s.size(), s) == 0;
    };
    const std::string w = at_line(source, ln);
    const auto value = [&](std::size_t i) { return to_number(token(t, i, w), w); };
    const auto whole = [&](std::size_t i) {
      const double v = value(i);
      if (v < 0.0 || v != std::floor(v)) throw std::invalid_argument(w + ": expected a count");
      return static_cast<std::size_t>(v);
    };
    if (ends("atom types")) {
      n_types = whole(0);
    } else if (ends("atoms")) {
      n_atoms = whole(0);
    } else if (ends("xy xz yz")) {
      xy = value(0);
      xz = value(1);
      yz = value(2);
    } else {
      static constexpr std::string_view box_keys[3] = {"xlo xhi", "ylo yhi", "zlo zhi"};
      for (int k = 0; k < 3; ++k) {
        if (!ends(box_keys[k])) continue;
        lo[k] = value(0);
        hi[k] = value(1);
        have[k] = true;
      }
    }
  }
  if (!have[0] || !have[1] || !have[2]) throw std::invalid_argument(source + ": the box (xlo xhi, ylo yhi, zlo zhi) is missing");
  cell.lattice = Lattice::from_vectors_A(
      Mat3::from_rows({hi[0] - lo[0], 0, 0}, {xy, hi[1] - lo[1], 0}, {xz, yz, hi[2] - lo[2]}));
  std::map<int, unsigned> type_z = opt.lammps_types;
  std::map<int, double> type_mass;
  std::map<int, unsigned> type_comment;
  struct Row {
    std::vector<std::string> t;
    std::size_t line;
  };
  std::vector<Row> atom_rows;
  std::string style = opt.lammps_atom_style;
  while (ln < lines.size()) {
    const std::string raw = lines[ln];
    const std::string l = strip_comment(raw);
    if (l.empty()) {
      ++ln;
      continue;
    }
    const std::string section = l;
    const std::string comment = raw.find('#') == std::string::npos ? "" : trim(raw.substr(raw.find('#') + 1));
    ++ln;
    std::vector<Row> rows;
    while (ln < lines.size()) {
      const std::string d = strip_comment(lines[ln]);
      if (d.empty()) {
        ++ln;
        if (!rows.empty()) break;
        continue;
      }
      if (!numeric_start(d)) break;
      rows.push_back({split_ws(d), ln});
      // Masses: a "# Fe" comment names the element.
      if (section == "Masses") {
        const std::string &full = lines[ln];
        if (const auto h = full.find('#'); h != std::string::npos) {
          if (const auto z = Element::find_z(trim(full.substr(h + 1)))) {
            type_comment[static_cast<int>(to_number(rows.back().t[0], at_line(source, ln)))] = *z;
          }
        }
      }
      ++ln;
    }
    if (section == "Masses") {
      for (const auto &r : rows) {
        if (r.t.size() < 2) throw std::invalid_argument(at_line(source, r.line) + ": expected type and mass");
        type_mass[static_cast<int>(to_number(r.t[0], at_line(source, r.line)))] = to_number(r.t[1], at_line(source, r.line));
      }
    } else if (section == "Atoms") {
      if (!comment.empty() && style.empty()) style = lower(split_ws(comment).front());
      atom_rows = std::move(rows);
    }
  }
  if (atom_rows.empty()) throw std::invalid_argument(source + ": no Atoms section");
  if (n_atoms != 0 && atom_rows.size() != n_atoms) {
    throw std::invalid_argument(source + ": the header says " + std::to_string(n_atoms) + " atoms, the Atoms section has " +
                                std::to_string(atom_rows.size()));
  }
  // Column of the type and of x for each style.
  int c_type = 1, c_x = 2;
  if (style.empty()) {
    const std::size_t n = atom_rows.front().t.size();
    if (n == 5 || n == 8) style = "atomic";
    else if (n == 7 || n == 10) style = "full";
    else {
      throw std::invalid_argument(source + ": cannot tell the atom style from " + std::to_string(n) +
                                  " columns; add '# atomic' (or charge, full, molecular) after 'Atoms' or give the style");
    }
  }
  if (style == "atomic") { c_type = 1; c_x = 2; }
  else if (style == "charge") { c_type = 1; c_x = 3; }
  else if (style == "full") { c_type = 2; c_x = 4; }
  else if (style == "molecular" || style == "bond" || style == "angle") { c_type = 2; c_x = 3; }
  else throw std::invalid_argument(source + ": atom style '" + style + "' is not supported (atomic, charge, full, molecular, bond, angle)");
  const auto element_of_type = [&](int type, std::size_t line) -> unsigned {
    if (const auto it = type_z.find(type); it != type_z.end()) return it->second;
    if (const auto it = type_comment.find(type); it != type_comment.end()) return it->second;
    const auto m = type_mass.find(type);
    if (m == type_mass.end()) {
      throw std::invalid_argument(at_line(source, line) + ": atom type " + std::to_string(type) +
                                  " has no mass and no element; give the type -> element map");
    }
    unsigned best = 0;
    double d_best = 1e9;
    for (const auto &e : Element::all()) {
      const double d = std::abs(in_units(e.natural_mass(), atomic_mass_unit) - m->second);
      if (d < d_best) {
        d_best = d;
        best = e.z();
      }
    }
    if (d_best > 0.1) {
      throw std::invalid_argument(source + ": atom type " + std::to_string(type) + " (mass " + fmt(m->second, 4) +
                                  " amu) matches no standard atomic weight within 0.1 amu; give the type -> element map");
    }
    return best;
  };
  std::vector<RawAtom> atoms;
  std::map<unsigned, int> per_element;
  for (std::size_t k = 0; k < atom_rows.size(); ++k) {
    const auto &r = atom_rows[k];
    const std::string w = at_line(source, r.line);
    if (static_cast<int>(r.t.size()) < c_x + 3) throw std::invalid_argument(w + ": too few columns for style " + style);
    const int type = static_cast<int>(to_number(r.t[static_cast<std::size_t>(c_type)], w + ": type"));
    if (n_types != 0 && (type < 1 || static_cast<std::size_t>(type) > n_types)) {
      throw std::invalid_argument(w + ": atom type " + std::to_string(type) + " outside 1.." + std::to_string(n_types));
    }
    const unsigned z = element_of_type(type, r.line);
    const Vec3 r_A{to_number(r.t[static_cast<std::size_t>(c_x)], w) - lo[0],
                   to_number(r.t[static_cast<std::size_t>(c_x + 1)], w) - lo[1],
                   to_number(r.t[static_cast<std::size_t>(c_x + 2)], w) - lo[2]};
    atoms.push_back({symbol(z) + std::to_string(++per_element[z]), z, 1.0, cell.lattice.to_fractional(r_A),
                     SiteKind::lattice, k});
  }
  cell.sites = merge_atoms(cell.lattice, atoms, opt.merge_distance_A, source);
  return cell;
}

// ---------------------------------------------------------------------------
// Extended XYZ
// ---------------------------------------------------------------------------

std::map<std::string, std::string> xyz_keys(const std::string &line, const std::string &where) {
  std::map<std::string, std::string> kv;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i >= line.size()) break;
    std::size_t k = i;
    while (k < line.size() && line[k] != '=' && !std::isspace(static_cast<unsigned char>(line[k]))) ++k;
    const std::string key = lower(line.substr(i, k - i));
    if (k >= line.size() || line[k] != '=') {
      kv[key] = "T"; // a bare key is a true flag
      i = k;
      continue;
    }
    ++k;
    std::string value;
    if (k < line.size() && (line[k] == '"' || line[k] == '\'')) {
      const char q = line[k];
      const std::size_t e = line.find(q, k + 1);
      if (e == std::string::npos) throw std::invalid_argument(where + ": unterminated quote in the comment line");
      value = line.substr(k + 1, e - k - 1);
      i = e + 1;
    } else {
      std::size_t e = k;
      while (e < line.size() && !std::isspace(static_cast<unsigned char>(line[e]))) ++e;
      value = line.substr(k, e - k);
      i = e;
    }
    kv[key] = value;
  }
  return kv;
}

UnitCell parse_extxyz(std::string_view text, const CellReadOptions &opt, const std::string &source) {
  const auto lines = lines_of(text);
  if (lines.size() < 2) throw std::invalid_argument(source + ": too short for an XYZ file");
  const auto count_line = split_ws(lines[0]);
  const double nd = to_number(token(count_line, 0, at_line(source, 0) + ": number of atoms"),
                              at_line(source, 0) + ": number of atoms");
  if (nd < 1 || nd != std::floor(nd)) throw std::invalid_argument(at_line(source, 0) + ": the number of atoms must be a positive integer");
  const auto n = static_cast<std::size_t>(nd);
  const auto kv = xyz_keys(lines[1], at_line(source, 1));
  const auto lat = kv.find("lattice");
  if (lat == kv.end()) throw std::invalid_argument(at_line(source, 1) + ": no Lattice=\"...\" (a plain XYZ file has no cell)");
  const auto lv = split_ws(lat->second);
  if (lv.size() != 9) throw std::invalid_argument(at_line(source, 1) + ": Lattice needs nine numbers");
  Mat3 m;
  for (int i = 0; i < 9; ++i) m[i / 3][i % 3] = to_number(lv[static_cast<std::size_t>(i)], at_line(source, 1) + ": Lattice");
  UnitCell cell;
  cell.lattice = Lattice::from_vectors_A(m);
  cell.name = kv.count("comment") ? kv.at("comment") : "";
  // Properties=name:type:count:...
  std::string props = kv.count("properties") ? kv.at("properties") : "species:S:1:pos:R:3";
  std::vector<std::string> parts;
  {
    std::string cur;
    for (char c : props) {
      if (c == ':') {
        parts.push_back(cur);
        cur.clear();
      } else {
        cur += c;
      }
    }
    parts.push_back(cur);
  }
  if (parts.size() % 3 != 0) throw std::invalid_argument(at_line(source, 1) + ": malformed Properties");
  int c_species = -1, c_pos = -1, c_z = -1, col = 0;
  for (std::size_t k = 0; k < parts.size(); k += 3) {
    const std::string name = lower(parts[k]);
    const int count = static_cast<int>(to_number(parts[k + 2], at_line(source, 1) + ": Properties"));
    if (name == "species") c_species = col;
    else if (name == "pos") c_pos = col;
    else if (name == "z" || name == "numbers") c_z = col;
    col += count;
  }
  if (c_pos < 0 || (c_species < 0 && c_z < 0)) {
    throw std::invalid_argument(at_line(source, 1) + ": Properties must include species (or Z) and pos");
  }
  if (lines.size() < n + 2) throw std::invalid_argument(source + ": fewer atom lines than the count says");
  std::vector<RawAtom> atoms;
  std::map<unsigned, int> per_element;
  for (std::size_t k = 0; k < n; ++k) {
    const std::string w = at_line(source, k + 2);
    const auto t = split_ws(lines[k + 2]);
    if (static_cast<int>(t.size()) < col) throw std::invalid_argument(w + ": expected " + std::to_string(col) + " columns");
    unsigned z = 0;
    if (c_species >= 0) {
      const auto e = Element::find_z(t[static_cast<std::size_t>(c_species)]);
      if (!e) throw std::invalid_argument(w + ": unknown element '" + t[static_cast<std::size_t>(c_species)] + "'");
      z = *e;
    } else {
      const double zd = to_number(t[static_cast<std::size_t>(c_z)], w + ": Z");
      if (zd < 1.0 || zd > Element::max_z || zd != std::floor(zd)) {
        throw std::invalid_argument(w + ": Z must be 1 to " + std::to_string(Element::max_z));
      }
      z = static_cast<unsigned>(zd);
    }
    const Vec3 r{to_number(t[static_cast<std::size_t>(c_pos)], w), to_number(t[static_cast<std::size_t>(c_pos + 1)], w),
                 to_number(t[static_cast<std::size_t>(c_pos + 2)], w)};
    atoms.push_back({symbol(z) + std::to_string(++per_element[z]), z, 1.0, cell.lattice.to_fractional(r),
                     SiteKind::lattice, k});
  }
  cell.sites = merge_atoms(cell.lattice, atoms, opt.merge_distance_A, source);
  return cell;
}

// ---------------------------------------------------------------------------
// Writers
// ---------------------------------------------------------------------------

std::string write_cif(const UnitCell &cell) {
  std::ostringstream o;
  std::string block;
  for (char c : cell.name.empty() ? std::string("amcset") : cell.name) {
    block += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
  }
  o << "# Written by AMCSET\n";
  o << "data_" << block << "\n";
  const auto deg = [](angle_quantity a) { return fmt(in_units(a, degree), 6); };
  o << "_cell_length_a " << fmt(in_units(cell.lattice.a(), angstrom), 8) << "\n";
  o << "_cell_length_b " << fmt(in_units(cell.lattice.b(), angstrom), 8) << "\n";
  o << "_cell_length_c " << fmt(in_units(cell.lattice.c(), angstrom), 8) << "\n";
  o << "_cell_angle_alpha " << deg(cell.lattice.alpha()) << "\n";
  o << "_cell_angle_beta " << deg(cell.lattice.beta()) << "\n";
  o << "_cell_angle_gamma " << deg(cell.lattice.gamma()) << "\n";
  o << "_symmetry_space_group_name_H-M 'P 1'\n";
  o << "loop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n";
  o << "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n"
       "_atom_site_fract_z\n_atom_site_occupancy\n_amcset_site_kind\n";
  for (const auto &s : cell.sites) {
    if (s.occupants.empty()) {
      o << "# empty " << to_string(s.kind) << " site " << s.label << " at " << fmt(s.position.x) << " "
        << fmt(s.position.y) << " " << fmt(s.position.z) << " (not written)\n";
      continue;
    }
    for (const auto &occ : s.occupants) {
      std::string label = s.label.empty() ? symbol(occ.z) : s.label;
      for (auto &c : label) {
        if (std::isspace(static_cast<unsigned char>(c))) c = '_';
      }
      o << label << " " << symbol(occ.z) << " " << fmt(s.position.x) << " " << fmt(s.position.y) << " "
        << fmt(s.position.z) << " " << fmt_general(occ.fraction, 8) << " " << to_string(s.kind) << "\n";
    }
  }
  return o.str();
}

std::string write_poscar(const UnitCell &cell) {
  const auto order = species_order(cell, "POSCAR");
  std::ostringstream o;
  std::string title = cell.name.empty() ? std::string("AMCSET cell") : cell.name;
  std::replace(title.begin(), title.end(), '\n', ' ');
  o << title << "\n1.0\n";
  const Mat3 &m = cell.lattice.vectors_A();
  for (int i = 0; i < 3; ++i) o << "  " << fmt(m[i][0]) << " " << fmt(m[i][1]) << " " << fmt(m[i][2]) << "\n";
  for (unsigned z : order) o << " " << symbol(z);
  o << "\n";
  for (unsigned z : order) {
    std::size_t n = 0;
    for (const auto &s : cell.sites) {
      if (!s.occupants.empty() && s.occupants.front().z == z) ++n;
    }
    o << " " << n;
  }
  o << "\nDirect\n";
  for (unsigned z : order) {
    for (const auto &s : cell.sites) {
      if (s.occupants.empty() || s.occupants.front().z != z) continue;
      o << "  " << fmt(s.position.x) << " " << fmt(s.position.y) << " " << fmt(s.position.z) << "\n";
    }
  }
  return o.str();
}

std::string write_extxyz(const UnitCell &cell) {
  species_order(cell, "Extended XYZ"); // checks occupancy
  std::size_t n = 0;
  for (const auto &s : cell.sites) {
    if (!s.occupants.empty()) ++n;
  }
  std::ostringstream o;
  const Mat3 &m = cell.lattice.vectors_A();
  o << n << "\nLattice=\"";
  for (int i = 0; i < 9; ++i) o << (i ? " " : "") << fmt(m[i / 3][i % 3]);
  o << "\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"";
  if (!cell.name.empty()) {
    std::string c = cell.name;
    std::replace(c.begin(), c.end(), '"', '\'');
    o << " comment=\"" << c << "\"";
  }
  o << "\n";
  for (const auto &s : cell.sites) {
    if (s.occupants.empty()) continue;
    const Vec3 r = cell.lattice.to_cartesian_A(s.position);
    o << symbol(s.occupants.front().z) << " " << fmt(r.x) << " " << fmt(r.y) << " " << fmt(r.z) << "\n";
  }
  return o.str();
}

std::string write_lammps_data(const UnitCell &cell) {
  const auto order = species_order(cell, "LAMMPS data");
  // LAMMPS's restricted triclinic box: a along x, b in the x-y plane.
  const Mat3 &m = cell.lattice.vectors_A();
  const double lx = norm(m[0]);
  const Vec3 ah = m[0] / lx;
  const double xy = dot(m[1], ah);
  const double ly = std::sqrt(std::max(0.0, dot(m[1], m[1]) - xy * xy));
  const double xz = dot(m[2], ah);
  const double yz = (dot(m[1], m[2]) - xy * xz) / ly;
  const double lz = std::sqrt(std::max(0.0, dot(m[2], m[2]) - xz * xz - yz * yz));
  const Mat3 box = Mat3::from_rows({lx, 0, 0}, {xy, ly, 0}, {xz, yz, lz});
  std::size_t n = 0;
  for (const auto &s : cell.sites) {
    if (!s.occupants.empty()) ++n;
  }
  std::ostringstream o;
  // The first line of a data file is its title (LAMMPS skips it).
  std::string title = cell.name.empty() ? std::string("AMCSET cell") : cell.name;
  std::replace(title.begin(), title.end(), '\n', ' ');
  o << title << "\n\n";
  o << n << " atoms\n" << order.size() << " atom types\n\n";
  o << "0.0 " << fmt(lx) << " xlo xhi\n0.0 " << fmt(ly) << " ylo yhi\n0.0 " << fmt(lz) << " zlo zhi\n";
  if (std::abs(xy) > 1e-12 || std::abs(xz) > 1e-12 || std::abs(yz) > 1e-12) {
    o << fmt(xy) << " " << fmt(xz) << " " << fmt(yz) << " xy xz yz\n";
  }
  o << "\nMasses\n\n";
  for (std::size_t t = 0; t < order.size(); ++t) {
    o << t + 1 << " " << fmt(in_units(Element::by_z(order[t]).natural_mass(), atomic_mass_unit), 6) << " # "
      << symbol(order[t]) << "\n";
  }
  o << "\nAtoms # atomic\n\n";
  std::size_t id = 0;
  for (const auto &s : cell.sites) {
    if (s.occupants.empty()) continue;
    const auto type = std::find(order.begin(), order.end(), s.occupants.front().z) - order.begin() + 1;
    const Vec3 r = box.transposed() * s.position;
    o << ++id << " " << type << " " << fmt(r.x) << " " << fmt(r.y) << " " << fmt(r.z) << "\n";
  }
  return o.str();
}

} // namespace

std::string to_string(CellFormat f) {
  switch (f) {
  case CellFormat::cif: return "cif";
  case CellFormat::poscar: return "poscar";
  case CellFormat::lammps_data: return "lammps-data";
  case CellFormat::extxyz: return "extxyz";
  }
  return "cif";
}

CellFormat parse_cell_format(std::string_view s) {
  const std::string l = lower(s);
  if (l == "cif") return CellFormat::cif;
  if (l == "poscar" || l == "vasp" || l == "contcar") return CellFormat::poscar;
  if (l == "lammps" || l == "lammps-data" || l == "lammps_data" || l == "data" || l == "lmp") return CellFormat::lammps_data;
  if (l == "xyz" || l == "extxyz" || l == "ext-xyz") return CellFormat::extxyz;
  throw std::invalid_argument("unknown cell format '" + std::string(s) + "' (cif, poscar, lammps-data, extxyz)");
}

std::optional<CellFormat> format_from_path(const std::filesystem::path &path) {
  const std::string ext = lower(path.extension().string());
  const std::string stem = lower(path.filename().string());
  if (ext == ".cif") return CellFormat::cif;
  if (ext == ".vasp" || ext == ".poscar" || stem.rfind("poscar", 0) == 0 || stem.rfind("contcar", 0) == 0) {
    return CellFormat::poscar;
  }
  if (ext == ".data" || ext == ".lmp" || ext == ".lammps" || stem.rfind("data.", 0) == 0) return CellFormat::lammps_data;
  if (ext == ".xyz" || ext == ".extxyz") return CellFormat::extxyz;
  return std::nullopt;
}

UnitCell parse_cell(std::string_view text, CellFormat format, const CellReadOptions &options, const std::string &source) {
  UnitCell cell;
  switch (format) {
  case CellFormat::cif: cell = parse_cif(text, options, source); break;
  case CellFormat::poscar: cell = parse_poscar(text, options, source); break;
  case CellFormat::lammps_data: cell = parse_lammps_data(text, options, source); break;
  case CellFormat::extxyz: cell = parse_extxyz(text, options, source); break;
  }
  cell.validate(source, 0.0); // positions, elements, occupancy (close sites were merged)
  return cell;
}

UnitCell read_cell(const std::filesystem::path &path, std::optional<CellFormat> format, const CellReadOptions &options) {
  const auto f = format ? format : format_from_path(path);
  if (!f) {
    throw std::invalid_argument(path.string() + ": cannot tell the format from the name; give it (cif, poscar, "
                                                "lammps-data, extxyz)");
  }
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error(path.string() + ": cannot open the file");
  std::ostringstream ss;
  ss << in.rdbuf();
  if (in.bad()) throw std::runtime_error(path.string() + ": read error");
  UnitCell cell = parse_cell(ss.str(), *f, options, path.filename().string());
  if (cell.name.empty()) cell.name = path.stem().string();
  return cell;
}

std::string write_cell(const UnitCell &cell, CellFormat format) {
  switch (format) {
  case CellFormat::cif: return write_cif(cell);
  case CellFormat::poscar: return write_poscar(cell);
  case CellFormat::lammps_data: return write_lammps_data(cell);
  case CellFormat::extxyz: return write_extxyz(cell);
  }
  return {};
}

void save_cell(const UnitCell &cell, const std::filesystem::path &path, std::optional<CellFormat> format) {
  const auto f = format ? format : format_from_path(path);
  if (!f) throw std::invalid_argument(path.string() + ": cannot tell the format from the name");
  const std::string text = write_cell(cell, *f);
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error(path.string() + ": cannot write the file");
  out << text;
  if (!out) throw std::runtime_error(path.string() + ": write error");
}

std::pair<Mat3, Vec3> parse_symmetry_operation(std::string_view text) {
  std::vector<std::string> parts;
  std::string cur;
  for (char c : text) {
    if (c == ',') {
      parts.push_back(cur);
      cur.clear();
    } else if (!std::isspace(static_cast<unsigned char>(c))) {
      cur += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
  }
  parts.push_back(cur);
  if (parts.size() != 3) throw std::invalid_argument("symmetry operation '" + std::string(text) + "': expected three parts");
  Mat3 m;
  Vec3 t;
  for (int r = 0; r < 3; ++r) {
    m[r] = Vec3{};
    const std::string &p = parts[static_cast<std::size_t>(r)];
    if (p.empty()) throw std::invalid_argument("symmetry operation '" + std::string(text) + "': empty part");
    std::size_t i = 0;
    while (i < p.size()) {
      double sign = 1.0;
      if (p[i] == '+' || p[i] == '-') {
        sign = p[i] == '-' ? -1.0 : 1.0;
        ++i;
      }
      if (i >= p.size()) throw std::invalid_argument("symmetry operation '" + std::string(text) + "': dangling sign");
      // Optional number (integer, decimal or fraction) ...
      double coef = 1.0;
      bool has_num = false;
      std::size_t j = i;
      while (j < p.size() && (std::isdigit(static_cast<unsigned char>(p[j])) || p[j] == '.')) ++j;
      if (j > i) {
        coef = std::stod(p.substr(i, j - i));
        has_num = true;
        i = j;
        if (i < p.size() && p[i] == '/') {
          std::size_t k = i + 1;
          while (k < p.size() && std::isdigit(static_cast<unsigned char>(p[k]))) ++k;
          if (k == i + 1) throw std::invalid_argument("symmetry operation '" + std::string(text) + "': bad fraction");
          coef /= std::stod(p.substr(i + 1, k - i - 1));
          i = k;
        }
        if (i < p.size() && p[i] == '*') ++i;
      }
      // ... followed by a variable, or standing alone (a translation).
      if (i < p.size() && (p[i] == 'x' || p[i] == 'y' || p[i] == 'z')) {
        m[r][p[i] - 'x'] += sign * coef;
        ++i;
      } else if (has_num) {
        t[r] += sign * coef;
      } else {
        throw std::invalid_argument("symmetry operation '" + std::string(text) + "': cannot read '" + p + "'");
      }
    }
  }
  return {m, t};
}

std::optional<unsigned> element_of_symbol(std::string_view text) {
  std::string letters;
  for (char c : text) {
    if (!std::isalpha(static_cast<unsigned char>(c))) break;
    letters += c;
    if (letters.size() == 2) break;
  }
  if (letters.empty()) return std::nullopt;
  if (letters.size() == 2) {
    std::string two{static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))),
                    static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1])))};
    // Only an element symbol, not a name ("Fe", not "Iron").
    if (const auto z = Element::find_z(two); z && Element::by_z(*z).symbol() == two) return z;
  }
  const std::string one{static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0])))};
  if (const auto z = Element::find_z(one); z && Element::by_z(*z).symbol() == one) return z;
  return std::nullopt;
}

} // namespace crystal
} // namespace amcset
