#include <QApplication>
#include <QPropertyAnimation>
#include <QTimer>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DExtras/QCuboidMesh>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DRender/QCamera>
#include <Qt3DRender/QPointLight>

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);

  // Create the main 3D window
  Qt3DExtras::Qt3DWindow view;

  // Root entity of the scene
  Qt3DCore::QEntity *rootEntity = new Qt3DCore::QEntity;

  // Cuboid entity
  Qt3DCore::QEntity *cuboidEntity = new Qt3DCore::QEntity(rootEntity);

  // Cuboid mesh component
  Qt3DExtras::QCuboidMesh *cuboidMesh = new Qt3DExtras::QCuboidMesh;
  cuboidMesh->setXExtent(1.0f); // Set cuboid dimensions
  cuboidMesh->setYExtent(1.0f); // Set cuboid dimensions
  cuboidMesh->setZExtent(1.0f); // Set cuboid dimensions

  // Transform component for the cuboid
  Qt3DCore::QTransform *cuboidTransform = new Qt3DCore::QTransform;
  cuboidTransform->setTranslation(
      QVector3D(0, 0, 0)); // Position the cuboid at the origin

  // Material component for the cuboid
  Qt3DExtras::QPhongMaterial *cuboidMaterial = new Qt3DExtras::QPhongMaterial;
  // cuboidMaterial->setDiffuse(QColor(QRgb(0x00FF00))); // Green color
  cuboidMaterial->setAmbient(QColor(QRgb(0xFF0000))); // Green color

  // Add components to the cuboid entity
  cuboidEntity->addComponent(cuboidMesh);
  cuboidEntity->addComponent(cuboidTransform);
  cuboidEntity->addComponent(cuboidMaterial);

  // Camera
  Qt3DRender::QCamera *camera = view.camera();
  camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 1000.0f);
  camera->setPosition(QVector3D(0, 0, 10));  // Initial camera position
  camera->setViewCenter(QVector3D(0, 0, 0)); // Look at the origin

  // Light
  Qt3DCore::QEntity *lightEntity = new Qt3DCore::QEntity(rootEntity);
  Qt3DRender::QPointLight *pointLight =
      new Qt3DRender::QPointLight(lightEntity);
  pointLight->setColor(Qt::white);
  pointLight->setIntensity(10.0f);
  lightEntity->addComponent(pointLight);

  // Set the root entity of the scene
  view.setRootEntity(rootEntity);

  // Animation to rotate the camera around the cube
  QPropertyAnimation *cameraAnimation =
      new QPropertyAnimation(camera, "position");
  cameraAnimation->setDuration(10000); // 10 seconds for a full rotation
  cameraAnimation->setLoopCount(-1);   // Infinite loop

  // Define the circular path for the camera
  QVector3D center(0, 0, 0); // Center of the cube
  float radius = 10.0f;      // Radius of the circular path
  cameraAnimation->setKeyValueAt(
      0, QVector3D(radius, 0, 0)); // Start at (radius, 0, 0)
  cameraAnimation->setKeyValueAt(
      0.25, QVector3D(10.0, radius, 10.0)); // Move to (0, radius, 0)
  cameraAnimation->setKeyValueAt(
      0.5, QVector3D(-radius, 0, 0)); // Move to (-radius, 0, 0)
  cameraAnimation->setKeyValueAt(
      0.75, QVector3D(0, -radius, 0)); // Move to (0, -radius, 0)
  cameraAnimation->setKeyValueAt(
      1, QVector3D(radius, 0, 0)); // Return to (radius, 0, 0)

  // Start the animation
  cameraAnimation->start();

  // Show the window
  view.show();

  return app.exec();
}

// #include <QApplication>
// #include <QGroupBox>
// #include <QHBoxLayout>
// #include <QLabel>
// #include <QMainWindow>
// #include <QPaintEvent>
// #include <QPainter>
// #include <QPushButton>
// #include <QSpinBox>
// #include <QVBoxLayout>
// #include <QWidget>
// #include <limits>

// struct Line {
//   int startX;
//   int startY;
//   int endX;
//   int endY;
// };

// // Custom widget that acts as our canvas
// class Canvas : public QWidget {
// public:
//   Canvas(QWidget *parent = nullptr) : QWidget(parent) {
//     // Set a minimum size for our canvas
//     setMinimumSize(1000, 1000);

//     // Initialize line coordinates
//     lines = std::vector<Line>(); // Prevent segfaults when mouse drawing on
//                                  // program opening.
//     shouldDrawLine = false;
//   }

//   // Set the coordinates for the line and trigger a repaint
//   void setLineCoordinates(int x1, int y1, int x2, int y2) {
//     lines.emplace_back(Line{x1, y1, x2, y2});
//     shouldDrawLine = true;
//     update(); // Request a repaint
//   }

// protected:
//   // Override the paint event to draw on the canvas
//   void paintEvent(QPaintEvent *event) override {
//     Q_UNUSED(event);
//     QPainter painter(this);
//     painter.setRenderHint(QPainter::Antialiasing, true);

//     // Fill the background with white
//     painter.fillRect(rect(), Qt::white);

//     // Only draw the line if we've been told to
//     if (shouldDrawLine) {
//       // Set up the pen (color, width, style)
//       QPen pen(Qt::blue, 3, Qt::SolidLine);
//       painter.setPen(pen);

//       // Draw the line using our stored coordinates
//       for (auto &line : lines) {
//         painter.drawLine(line.startX, line.startY, line.endX, line.endY);
//       }
//     }
//   }

//   void mousePressEvent(QMouseEvent *event) override {
//     mousePos = event->pos();
//     const auto lastPos = lines.back();
//     lines.emplace_back(
//         Line{lastPos.endX, lastPos.endY, mousePos.x(), mousePos.y()});
//     update();
//   }

// private:
//   std::vector<Line> lines;
//   bool shouldDrawLine;
//   QPoint mousePos;
// };

// // Main window to contain our canvas and controls
// class MainWindow : public QMainWindow {
// public:
//   MainWindow() {
//     // Set window title
//     setWindowTitle("Qt Canvas Line Example");

//     // Create the main widget and layout
//     QWidget *centralWidget = new QWidget(this);
//     QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);

//     // Create canvas
//     canvas = new Canvas(this);

//     // Create input controls
//     QGroupBox *inputGroup = new QGroupBox("Line Coordinates", this);
//     QGridLayout *inputLayout = new QGridLayout(inputGroup);

//     // Start point controls
//     inputLayout->addWidget(new QLabel("Start X:"), 0, 0);
//     startXSpinBox = new QSpinBox(inputGroup);
//     startXSpinBox->setRange(0, std::numeric_limits<int>::max());
//     inputLayout->addWidget(startXSpinBox, 0, 1);

//     inputLayout->addWidget(new QLabel("Start Y:"), 0, 2);
//     startYSpinBox = new QSpinBox(inputGroup);
//     startYSpinBox->setRange(0, std::numeric_limits<int>::max());
//     inputLayout->addWidget(startYSpinBox, 0, 3);

//     // End point controls
//     inputLayout->addWidget(new QLabel("End X:"), 1, 0);
//     endXSpinBox = new QSpinBox(inputGroup);
//     endXSpinBox->setRange(0, std::numeric_limits<int>::max());
//     inputLayout->addWidget(endXSpinBox, 1, 1);

//     inputLayout->addWidget(new QLabel("End Y:"), 1, 2);
//     endYSpinBox = new QSpinBox(inputGroup);
//     endYSpinBox->setRange(0, std::numeric_limits<int>::max());
//     inputLayout->addWidget(endYSpinBox, 1, 3);

//     // Draw button
//     QPushButton *drawButton = new QPushButton("Draw Line", inputGroup);
//     inputLayout->addWidget(drawButton, 2, 0, 1, 4);

//     // Connect the button's clicked signal to our slot
//     connect(drawButton, &QPushButton::clicked, this,
//             &MainWindow::onDrawButtonClicked);

//     // Add the widgets to the main layout
//     mainLayout->addWidget(canvas);
//     mainLayout->addWidget(inputGroup);

//     // Set the central widget
//     setCentralWidget(centralWidget);

//     onDrawButtonClicked(); // Add a dummy coordinate to prevent empty vector
//                            // segfaults
//   }

// private slots:
//   void onDrawButtonClicked() {
//     // Get the coordinates from the spin boxes
//     int x1 = startXSpinBox->value();
//     int y1 = startYSpinBox->value();
//     int x2 = endXSpinBox->value();
//     int y2 = endYSpinBox->value();

//     // Tell the canvas to draw the line
//     canvas->setLineCoordinates(x1, y1, x2, y2);
//   }

// private:
//   Canvas *canvas;
//   QSpinBox *startXSpinBox, *startYSpinBox, *endXSpinBox, *endYSpinBox;
// };

// int main(int argc, char *argv[]) {
//   // Initialize the Qt application
//   QApplication app(argc, argv);

//   // Create and show the main window
//   MainWindow window;
//   window.show();

//   // Enter the application event loop
//   return app.exec();
// }

// #include <atomic>
// #include <boost/asio.hpp>
// #include <boost/asio/executor_work_guard.hpp>
// #include <boost/lockfree/spsc_queue.hpp>
// #include <chrono>
// #include <cstdlib>
// #include <fstream>
// #include <iostream>
// #include <memory>
// #include <random>
// #include <thread>
// #include <variant>
// #include <vector>

// using DataVariant = std::variant<std::vector<double>, int, double>;
// using namespace boost::asio;

// struct ThreadTerminator {};
// // using QueueVariant =
// //     std::variant<DataVariant, ThreadTerminator>; // This can be simplified
// to
// //     be
// // a single variant

// boost::lockfree::spsc_queue<DataVariant *> queue(1024 *
//                                                  10); // Set size of 1024
//                                                  bytes
// std::atomic<bool> should_terminate{false};            // Atomic access bool
// std::atomic<size_t> pending_items{0};

// struct VariantVisitor {
//   VariantVisitor(const std::string &filename)
//       : output_file(filename, std::ios::app) {
//     std::cout << "Opening file: " << filename << std::endl;
//     if (!output_file.is_open()) {
//       throw std::runtime_error("Failed to open output file: " + filename);
//     }
//   }

//   void operator()(const DataVariant &data) const {
//     std::cout << "Processing DataVariant" << std::endl;
//     std::visit([this](const auto &value) { this->process_value(value); },
//     data);
//   }

//   // void operator()(const ThreadTerminator &) const {
//   //   std::cout << "Processing ThreadTerminator" << std::endl;
//   // }

//   ~VariantVisitor() {
//     std::cout << "Closing file" << std::endl;
//     if (output_file.is_open()) {
//       output_file.close();
//     }
//   }

// private:
//   void process_value(const std::vector<double> &vec) const {
//     std::cout << "Writing vector to file" << std::endl;
//     output_file << "Vector: ";
//     for (auto v : vec) {
//       output_file << v << " ";
//     }
//     output_file << std::endl;
//     output_file.flush();
//   }

//   void process_value(int x) const {
//     std::cout << "Writing int to file" << std::endl;
//     output_file << "Int: " << x << std::endl;
//     output_file.flush();
//   }

//   void process_value(double x) const {
//     std::cout << "Writing double to file" << std::endl;
//     output_file << "Double: " << x << std::endl;
//     output_file.flush();
//   }

//   mutable std::ofstream output_file;
// };

// void output_thread_func() {
//   std::cout << "Output thread started" << std::endl;
//   VariantVisitor visitor("/home/nathaniel/Code/C++/AMCSET/test_output.txt");

//   while (!should_terminate || pending_items > 0) {
//     DataVariant *raw_data = nullptr;
//     if (queue.pop(raw_data)) {
//       std::cout << "Popped item from queue, pending_items=" << pending_items
//                 << std::endl;
//       std::unique_ptr<DataVariant> data(raw_data);

//       // if (!std::holds_alternative<ThreadTerminator>(*data)) {
//       std::cout << "Processing data item" << std::endl;
//       std::visit(visitor, *data);
//       pending_items--;
//       std::cout << "Processed item, remaining pending_items=" <<
//       pending_items
//                 << std::endl;
//       // free(raw_data);
//       // }
//     } else {
//       std::this_thread::yield();
//     }
//   }
//   std::cout << "Output thread terminating" << std::endl;
// }

// auto rng = std::mt19937();
// auto dist = std::uniform_real_distribution<double>();
// auto rng_mutex = std::mutex();

// DataVariant *bombardment_task() {
//   std::cout << "Task starting" << std::endl;
//   std::this_thread::sleep_for(std::chrono::milliseconds(100));
//   std::cout << "Task finished" << std::endl;
//   double test;
//   while (true) {
//     if (rng_mutex.try_lock()) {
//       test = dist(rng);
//       break;
//     }
//   }
//   rng_mutex.unlock();

//   if (test < 0.3) {
//     return new DataVariant((int)(test * 300.0));
//   } else if (test < 0.6) {
//     return new DataVariant(test);
//   } else {
//     std::vector<double> res;
//     res.push_back(dist(rng));
//     res.push_back(dist(rng));
//     res.push_back(dist(rng));
//     return new DataVariant(res);
//   }
// }

// void simulation(io_context &io) {
//   std::vector<std::thread> pool; // Spawn a threadpool

//   const int num_threads = std::thread::hardware_concurrency();
//   const int tasks_per_thread = 10;

//   // Use a work guard to keep the io_context active
//   auto work = boost::asio::make_work_guard(io);

//   std::cout << "Starting simulation threads" << std::endl;
//   for (int i = 0; i < num_threads; ++i) {
//     pool.emplace_back([&io]() {
//       for (int j = 0; j < tasks_per_thread; ++j) {
//         std::cout << "Thread " << std::this_thread::get_id()
//                   << " starting task " << j << std::endl;
//         auto result = bombardment_task();

//         post(io, [result]() {
//           std::cout << "Post callback executing, incrementing pending_items"
//                     << std::endl;
//           pending_items++;
//           bool success = queue.push(result);
//           if (!success) {
//             std::cerr << "Failed to push to queue!" << std::endl;
//             pending_items--;
//           }
//           std::cout << "Pushed to queue, pending_items=" << pending_items
//                     << std::endl;
//         });
//       }
//     });
//   }

//   // Wait for all pool threads to finish posting tasks
//   std::cout << "Joining pool threads" << std::endl;
//   for (auto &thread : pool) {
//     if (thread.joinable()) {
//       thread.join();
//     }
//   }

//   // Reset the work guard to allow the io_context to exit
//   work.reset();

//   std::cout << "Running io_context" << std::endl;
//   io.run();

//   std::cout << "All threads joined, sending terminator" << std::endl;
//   should_terminate = true;
//   // auto *terminator = new QueueVariant(ThreadTerminator{});
//   // queue.push(terminator);
// }

// int main() {
//   std::cout << "Program starting" << std::endl;
//   io_context io;

//   std::thread output_thread(output_thread_func);

//   simulation(io);

//   std::cout << "Waiting for output thread" << std::endl;
//   if (output_thread.joinable()) {
//     output_thread.join();
//   }

//   std::cout << "Program completed" << std::endl;
//   return 0;
// }
