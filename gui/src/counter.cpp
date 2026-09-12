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
