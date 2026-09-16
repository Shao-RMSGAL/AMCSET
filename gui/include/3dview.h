#pragma once

#include <QMutex>
#include <QVBoxLayout>
#include <QVector3D>
#include <QVector>
#include <QWidget>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DExtras/QCuboidMesh>
#include <Qt3DExtras/QOrbitCameraController>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QSphereMesh>
#include <Qt3DExtras/Qt3DWindow>
// #include <Qt3DRender/QAttribute>
// #include <Qt3DRender/QBuffer>
#include <Qt3DRender/QCamera>
// #include <Qt3DRender/QGeometry>
#include <QQueue>
#include <Qt3DRender/QGeometryRenderer>
#include <Qt3DRender/QMesh>

class Widget3DView : public QWidget {
  Q_OBJECT

public:
  explicit Widget3DView(QWidget *parent = nullptr);
  ~Widget3DView();

  // Function to draw a line between two 3D points
  void drawLine(const QVector3D &from, const QVector3D &to,
                const QColor &color = Qt::red);

  // Get root entity - allows adding additional entities if needed
  Qt3DCore::QEntity *rootEntity() const { return m_rootEntity; }

  // Get the 3D view - allows adjusting camera or other view settings
  Qt3DExtras::Qt3DWindow *view() const { return m_view; }

  void requestDrawLine(const QVector3D &from, const QVector3D &to,
                       const QColor &color = Qt::red);

  void clearScene();
  void requestClearScene();
  void setRange(double range);

private slots:
  // Process pending line drawing requests
  void processPendingLineRequests();

private:
  Qt3DExtras::Qt3DWindow *m_view;
  QWidget *m_container;
  Qt3DCore::QEntity *m_rootEntity;

  // Helper method to create a line entity
  Qt3DCore::QEntity *createLineEntity(const QVector3D &from,
                                      const QVector3D &to, const QColor &color);
  struct LineRequest {
    QVector3D from;
    QVector3D to;
    QColor color;
  };

  bool m_clearSceneRequested = false;
  QQueue<LineRequest> m_lineRequestQueue;
  QMutex m_lineRequestMutex;
  QTimer *m_updateTimer;
  QList<Qt3DCore::QEntity *> m_sceneEntities;
};
