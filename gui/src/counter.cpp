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

#include "counter.h"

Counter::Counter(QObject *parent)
    : QObject(parent), m_value(0), m_changed(false) {}

int Counter::getValue() const {
  return m_value.load(std::memory_order_acquire);
}

void Counter::setValue(int value) {
  const auto previous = m_value.exchange(value, std::memory_order_acq_rel);
  if (previous != value) {
    m_changed.store(true, std::memory_order_release);
    emit valueChanged(value);
  }
}

void Counter::increment() {
  const auto next = m_value.fetch_add(1, std::memory_order_acq_rel) + 1;
  m_changed.store(true, std::memory_order_release);
  emit valueChanged(next);
}

void Counter::reset() {
  m_value.store(0, std::memory_order_release);
  m_changed.store(false, std::memory_order_release);
  emit valueChanged(0);
}

bool Counter::hasChanged() {
  return m_changed.exchange(false, std::memory_order_acq_rel);
}
