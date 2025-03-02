#pragma once
#include <QMutex>
#include <QObject>

#include <boost/thread/mutex.hpp>

class Counter : public QObject {
  Q_OBJECT
public:
  explicit Counter(QObject *parent = nullptr);

  int getValue() const;
  void increment();
  bool hasChanged();
  void reset();

signals:
  void valueChanged(int newValue);

private:
  mutable boost::mutex m_mutex;
  mutable QMutex q_mutex;
  int m_value;
  mutable bool m_changed;
};
