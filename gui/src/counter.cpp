#include "counter.h"

// FIX: Does not work. Abandon for queue approach
Counter::Counter(QObject *parent)
    : QObject(parent), m_value(0), m_changed(false) {}

int Counter::getValue() const {
  QMutexLocker locker(&q_mutex);
  m_mutex.lock();
  m_changed = false;
  m_mutex.unlock();
  return m_value;
}

void Counter::increment() {
  QMutexLocker locker(&q_mutex);
  m_mutex.lock();
  m_value++;
  m_changed = true;
  m_mutex.unlock();
  emit valueChanged(m_value);
}

void Counter::reset() {
  QMutexLocker locker(&q_mutex);
  m_mutex.lock();
  m_value = 0;
  m_changed = false;
  m_mutex.unlock();
  emit valueChanged(m_value);
}

bool Counter::hasChanged() {
  QMutexLocker locker(&q_mutex);
  m_mutex.lock();
  auto result = m_changed;
  m_mutex.unlock();
  return result;
}
