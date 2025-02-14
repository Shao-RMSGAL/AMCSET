#include <QWidget>

class Canvas : public QWidget {
  Q_OBJECT

public:
  Canvas(QWidget *parent = nullptr);

protected:
  void paintEvent(QPaintEvent *event) override;
};
