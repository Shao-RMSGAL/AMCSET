#include "canvas.h"
#include <QPainter>
#include <QPainterPath>

Canvas::Canvas(QWidget *parent) : QWidget(parent) {
  // Set white background
  setAutoFillBackground(true);
  QPalette pal = palette();
  pal.setColor(QPalette::Window, Qt::white);
  setPalette(pal);
}

void Canvas::paintEvent(QPaintEvent *) {
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);

  // Draw a house

  // Main house structure (rectangle)
  painter.setPen(QPen(Qt::black, 2));
  painter.setBrush(QBrush(QColor(255, 220, 180)));
  painter.drawRect(300, 200, 200, 150);

  // Roof (triangle)
  QPolygon roof;
  roof << QPoint(280, 200) << QPoint(400, 120) << QPoint(520, 200);
  painter.setBrush(QBrush(QColor(139, 69, 19)));
  painter.drawPolygon(roof);

  // Door
  painter.setBrush(QBrush(QColor(139, 69, 19)));
  painter.drawRect(375, 270, 50, 80);

  // Doorknob
  painter.setBrush(QBrush(Qt::yellow));
  painter.drawEllipse(410, 310, 8, 8);

  // Windows
  painter.setBrush(QBrush(Qt::lightGray));
  painter.drawRect(320, 230, 40, 40); // Left window
  painter.drawRect(440, 230, 40, 40); // Right window

  // Window crossbars
  painter.drawLine(320, 250, 360, 250); // Left window horizontal
  painter.drawLine(340, 230, 340, 270); // Left window vertical
  painter.drawLine(440, 250, 480, 250); // Right window horizontal
  painter.drawLine(460, 230, 460, 270); // Right window vertical

  // Sun
  painter.setPen(Qt::NoPen);
  painter.setBrush(QBrush(Qt::yellow));
  painter.drawEllipse(50, 50, 60, 60);

  // Ground
  painter.setPen(QPen(Qt::black, 2));
  painter.setBrush(QBrush(QColor(34, 139, 34)));
  painter.drawRect(0, 350, width(), 250);
}
