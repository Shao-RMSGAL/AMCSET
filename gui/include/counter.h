#pragma once
#include <QMutex>
#include <QObject>

#include <atomic>
#include <boost/thread/mutex.hpp>

class Counter : public QObject {
  Q_OBJECT
public:
  explicit Counter(QObject *parent = nullptr);

  int getValue() const;
  void increment();
  void setValue(int value);
  bool hasChanged();
  void reset();

signals:
  void valueChanged(int newValue);

private:
  std::atomic<int> m_value{0};
  std::atomic<bool> m_changed{false};
};
