#include "3dview.h"

#include <QTimer>

Widget3DView::Widget3DView(QWidget *parent) : QWidget(parent) {
  // Create layout for widget
  QVBoxLayout *layout = new QVBoxLayout(this);
  layout->setContentsMargins(0, 0, 0, 0);

  // Create 3D window
  m_view = new Qt3DExtras::Qt3DWindow();
  // m_view->defaultFrameGraph()->setClearColor(QColor(50, 50, 50));

  // Create container widget for the 3D window
  m_container = QWidget::createWindowContainer(m_view, this);
  m_container->setMinimumSize(200, 200);
  m_container->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

  // Set up the layout with our container widget
  layout->addWidget(m_container);

  // Root entity
  m_rootEntity = new Qt3DCore::QEntity();

  // Set up camera
  auto *camera = m_view->camera();
  auto defaultRange = 1000.0f; // TODO: Use a single source of truth for this
  camera->setPosition(QVector3D(defaultRange, 0.0f, -defaultRange));
  camera->setViewCenter(QVector3D(0, 0, 0));
  // Increase the far clipping plane distance
  auto *lens = camera->lens();
  lens->setFarPlane(defaultRange * 10); // Adjust value as needed

  // You may also want to adjust the near plane if needed
  lens->setNearPlane(0.1f);

  // Camera controller for user interaction
  Qt3DExtras::QOrbitCameraController *camController =
      new Qt3DExtras::QOrbitCameraController(m_rootEntity);
  camController->setCamera(camera);

  // Set the root entity to the view
  m_view->setRootEntity(m_rootEntity);
  m_updateTimer = new QTimer(this);
  m_updateTimer->setInterval(16); // ~60 FPS
  connect(m_updateTimer, &QTimer::timeout, this,
          &Widget3DView::processPendingLineRequests);
  m_updateTimer->start();
}

void Widget3DView::clearScene() {
  // Delete all entities we've added to the scene
  foreach (Qt3DCore::QEntity *entity, m_sceneEntities) {
    // Remove from parent first (ensures proper cleanup in Qt3D)
    entity->setParent((Qt3DCore::QNode *)nullptr);
    entity->deleteLater();
  }

  // Clear the list
  m_sceneEntities.clear();
}

void Widget3DView::setRange(double range) {
  auto *camera = m_view->camera();
  camera->setPosition(QVector3D(range, 0.0f, -range));
  camera->lens()->setFarPlane(range * 10); // Adjust value as needed
}

void Widget3DView::requestDrawLine(const QVector3D &from, const QVector3D &to,
                                   const QColor &color) {
  // This method can be called from any thread
  // m_mutex.lock();
  m_lineRequestQueue.enqueue({from, to, color});
  // m_mutex.unlock();
}
void Widget3DView::requestClearScene() {
  // This method can be called from any thread
  // m_mutex.lock();
  m_clearSceneRequested = true;
  // m_mutex.unlock();
}

void Widget3DView::processPendingLineRequests() {
  // This method runs in the UI thread
  // m_mutex.lock();

  while (!m_lineRequestQueue.isEmpty()) {
    LineRequest request = m_lineRequestQueue.dequeue();
    // m_mutex.unlock(); // Unlock while drawing to prevent blocking other
    // threads

    // Draw the line (this must happen in the UI thread)
    drawLine(request.from, request.to, request.color);

    // m_mutex.lock();
  }

  // m_mutex.unlock();
}

Widget3DView::~Widget3DView() {
  if (m_updateTimer) {
    m_updateTimer->stop();
    m_updateTimer->deleteLater();
    m_updateTimer = nullptr;
  }

  clearScene();

  if (m_rootEntity) {
    m_rootEntity->deleteLater();
    m_rootEntity = nullptr;
  }

  if (m_view) {
    m_view->deleteLater();
    m_view = nullptr;
  }

  if (m_container) {
    m_container->deleteLater();
    m_container = nullptr;
  }
}

void Widget3DView::drawLine(const QVector3D &from, const QVector3D &to,
                            const QColor &color) {
  // Create a line entity and add it to the root entity
  Qt3DCore::QEntity *lineEntity = createLineEntity(from, to, color);
  lineEntity->setParent(m_rootEntity);
  m_sceneEntities.append(lineEntity);
}

Qt3DCore::QEntity *Widget3DView::createLineEntity(const QVector3D &from,
                                                  const QVector3D &to,
                                                  const QColor &color) {
  // Create entity for the line
  Qt3DCore::QEntity *lineEntity = new Qt3DCore::QEntity();

  // Create geometry
  Qt3DCore::QGeometry *geometry = new Qt3DCore::QGeometry(lineEntity);

  // Create position buffer
  QVector<float> vertices;
  vertices << from.x() << from.y() << from.z();
  vertices << to.x() << to.y() << to.z();

  Qt3DCore::QBuffer *vertexBuffer = new Qt3DCore::QBuffer(geometry);
  vertexBuffer->setData(
      QByteArray(reinterpret_cast<const char *>(vertices.constData()),
                 vertices.size() * sizeof(float)));

  // Create position attribute
  Qt3DCore::QAttribute *positionAttribute = new Qt3DCore::QAttribute(geometry);
  positionAttribute->setName(
      Qt3DCore::QAttribute::defaultPositionAttributeName());
  positionAttribute->setVertexBaseType(Qt3DCore::QAttribute::Float);
  positionAttribute->setVertexSize(3);
  positionAttribute->setAttributeType(Qt3DCore::QAttribute::VertexAttribute);
  positionAttribute->setBuffer(vertexBuffer);
  positionAttribute->setByteStride(3 * sizeof(float));
  positionAttribute->setCount(2);
  geometry->addAttribute(positionAttribute);

  // Create index buffer
  QVector<unsigned int> indices = {0, 1};
  Qt3DCore::QBuffer *indexBuffer = new Qt3DCore::QBuffer(geometry);

  indexBuffer->setData(
      QByteArray(reinterpret_cast<const char *>(indices.constData()),
                 indices.size() * sizeof(unsigned int)));

  // Create index attribute
  Qt3DCore::QAttribute *indexAttribute = new Qt3DCore::QAttribute(geometry);
  indexAttribute->setVertexBaseType(Qt3DCore::QAttribute::UnsignedInt);
  indexAttribute->setAttributeType(Qt3DCore::QAttribute::IndexAttribute);
  indexAttribute->setBuffer(indexBuffer);
  indexAttribute->setCount(2);
  geometry->addAttribute(indexAttribute);

  // Create geometry renderer
  Qt3DRender::QGeometryRenderer *renderer =
      new Qt3DRender::QGeometryRenderer(lineEntity);
  renderer->setGeometry(geometry);
  renderer->setPrimitiveType(Qt3DRender::QGeometryRenderer::Lines);
  lineEntity->addComponent(renderer);

  // Create material
  Qt3DExtras::QPhongMaterial *material =
      new Qt3DExtras::QPhongMaterial(lineEntity);
  material->setAmbient(color);
  lineEntity->addComponent(material);

  return lineEntity;
}
