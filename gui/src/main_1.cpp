
//=========================SPSC================================================
#include <atomic>
#include <boost/asio.hpp>
#include <boost/asio/executor_work_guard.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <thread>
#include <variant>
#include <vector>

using DataVariant = std::variant<std::vector<double>, int, double>;
using namespace boost::asio;

constexpr size_t num_bomb = 100;
boost::lockfree::spsc_queue<DataVariant *> queue(sizeof(DataVariant) *
                                                 num_bomb);
std::atomic<bool> should_terminate{false};
std::atomic<size_t> pending_items{0};

struct VariantVisitor {
  VariantVisitor(const std::string &filename)
      : output_file(filename, std::ios::app) {
    std::cout << "Opening file: " << filename << std::endl;
    if (!output_file.is_open()) {
      throw std::runtime_error("Failed to open output file: " + filename);
    }
  }

  void operator()(const DataVariant &data) const {
    std::cout << "Processing DataVariant" << std::endl;
    std::visit([this](const auto &value) { this->process_value(value); }, data);
  }

  ~VariantVisitor() {
    std::cout << "Closing file" << std::endl;
    if (output_file.is_open()) {
      output_file.close();
    }
  }

private:
  void process_value(const std::vector<double> &vec) const {
    std::cout << "Writing vector to file" << std::endl;
    output_file << "Vector: ";
    for (auto v : vec) {
      output_file << v << " ";
    }
    output_file << std::endl;
    output_file.flush();
  }

  void process_value(int x) const {
    std::cout << "Writing int to file" << std::endl;
    output_file << "Int: " << x << std::endl;
    output_file.flush();
  }

  void process_value(double x) const {
    std::cout << "Writing double to file" << std::endl;
    output_file << "Double: " << x << std::endl;
    output_file.flush();
  }

  mutable std::ofstream output_file;
};

void output_thread_func() {
  std::cout << "Output thread started" << std::endl;
  VariantVisitor visitor("/home/nathaniel/Code/C++/AMCSET/test_output.txt");

  while (!should_terminate || pending_items > 0) {
    DataVariant *raw_data = nullptr;
    if (queue.pop(raw_data)) {
      std::cout << "Popped item from queue, pending_items=" << pending_items
                << std::endl;
      std::unique_ptr<DataVariant> data(raw_data);
      std::cout << "Processing data item" << std::endl;
      std::visit(visitor, *data);
      pending_items--;
      std::cout << "Processed item, remaining pending_items=" << pending_items
                << std::endl;
    } else {
      std::this_thread::yield();
    }
  }
  std::cout << "Output thread terminating" << std::endl;
}

auto rng = std::mt19937();
auto dist = std::uniform_real_distribution<double>();
auto rng_mutex = std::mutex();

DataVariant *bombardment_task() {
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  double test;
  double attempt = 0;
  while (true) {
    attempt++;
    if (rng_mutex.try_lock()) {
      test = dist(rng);
      break;
    }
  }
  rng_mutex.unlock();

  if (test < 0.3) {
    return new DataVariant((int)(test * 300.0));
  } else if (test < 0.6) {
    return new DataVariant(test);
  } else {
    std::vector<double> res;
    attempt = 0;
    while (true) {
      attempt++;
      if (rng_mutex.try_lock()) {
        res.push_back(dist(rng));
        res.push_back(dist(rng));
        res.push_back(dist(rng));
        break;
      }
    }
    rng_mutex.unlock();
    return new DataVariant(res);
  }
}

std::atomic<size_t> threads_done = 0;
void simulation(io_context &io) {
  std::vector<std::thread> pool; // Spawn a threadpool

  const int num_threads = 14; // std::thread::hardware_concurrency();
  const int tasks_per_thread = 10;

  // Use a work guard to keep the io_context active
  auto work = boost::asio::make_work_guard(io);

  std::cout << "Starting simulation threads" << std::endl;
  for (int i = 0; i < num_threads; ++i) {
    pool.emplace_back([&io]() {
      for (int j = 0; j < tasks_per_thread; ++j) {
        std::cout << "Thread " << std::this_thread::get_id()
                  << " starting task " << j << std::endl;
        auto result = bombardment_task();

        post(io, [result]() {
          std::cout << "Post callback executing, incrementing pending_items"
                    << std::endl;
          pending_items++;
          bool success = queue.push(result);
          if (!success) {
            std::cerr << "Failed to push to queue!" << std::endl;
            pending_items--;
          } // If push fails, retry, don't throw away
          std::cout << "Pushed to queue, pending_items=" << pending_items
                    << std::endl;
          threads_done++;
          std::cout << "Threads done: " << threads_done << std::endl;
        });
        std::cout << "Bombardment thread " << std::this_thread::get_id()
                  << " exiting" << std::endl;
      }
    });
  }

  std::cout << "Starting IO Context" << std::endl;
  // work.reset();
  std::thread io_thread([&io]() { io.run(); });
  std::cout << "Finished IO Context" << std::endl;

  // Wait for all pool threads to finish posting tasks
  std::cout << "Joining pool threads" << std::endl;
  auto num_joined = 0;
  for (auto &thread : pool) {
    std::cout << "Joinable?" << std::endl;
    if (thread.joinable()) {
      std::cout << "Trying join" << std::endl;
      thread.join();
      num_joined++;
      std::cout << "Joined " << num_joined << std::endl;
    } else {
      std::cout << "Not Joinable." << std::endl;
    }
  }

  std::cout << "Pool threads joined" << std::endl;

  should_terminate = true;

  // Reset the work guard to allow the io_context to exit

  std::cout << "Running io_context" << std::endl;

  std::cout << "All threads joined, sending terminator" << std::endl;

  while (true) { // HACK:: Spin lock, do better
    std::cout << "Spinning: " << threads_done << std::endl;
    if (threads_done == 140) {
      io.stop();
      break;
    }
  }

  if (io_thread.joinable())
    io_thread.join();
}

int main() {
  std::cout << "Program starting..." << std::endl;
  io_context io;

  std::thread output_thread(output_thread_func);
  simulation(io);

  std::cout << "Waiting for output thread" << std::endl;
  if (output_thread.joinable()) {
    output_thread.join();
  }

  std::cout << "Program completed. Size of std::vector<double>: "
            << sizeof(std::vector<double>)
            << ". Variant: " << sizeof(DataVariant) << std::endl;
  return 0;
}
//=======================END
// SPSQ==============================================

// //===========================LINE
// DRAWING====================================== #include <QApplication>
// #include <QPropertyAnimation>
// #include <QTimer>
// #include <Qt3DCore/QEntity>
// #include <Qt3DCore/QTransform>
// #include <Qt3DExtras/QCuboidMesh>
// #include <Qt3DExtras/QPhongMaterial>
// #include <Qt3DExtras/Qt3DWindow>
// #include <Qt3DRender/QCamera>
// #include <Qt3DRender/QPointLight>

// int main(int argc, char *argv[]) {
//   QApplication app(argc, argv);

//   // Create the main 3D window
//   Qt3DExtras::Qt3DWindow view;

//   // Root entity of the scene
//   Qt3DCore::QEntity *rootEntity = new Qt3DCore::QEntity;

//   // Cuboid entity
//   Qt3DCore::QEntity *cuboidEntity = new Qt3DCore::QEntity(rootEntity);

//   // Cuboid mesh component
//   Qt3DExtras::QCuboidMesh *cuboidMesh = new Qt3DExtras::QCuboidMesh;
//   cuboidMesh->setXExtent(1.0f); // Set cuboid dimensions
//   cuboidMesh->setYExtent(1.0f); // Set cuboid dimensions
//   cuboidMesh->setZExtent(1.0f); // Set cuboid dimensions

//   // Transform component for the cuboid
//   Qt3DCore::QTransform *cuboidTransform = new Qt3DCore::QTransform;
//   cuboidTransform->setTranslation(
//       QVector3D(0, 0, 0)); // Position the cuboid at the origin

//   // Material component for the cuboid
//   Qt3DExtras::QPhongMaterial *cuboidMaterial = new
//   Qt3DExtras::QPhongMaterial;
//   // cuboidMaterial->setDiffuse(QColor(QRgb(0x00FF00))); // Green color
//   cuboidMaterial->setAmbient(QColor(QRgb(0xFF0000))); // Green color

//   // Add components to the cuboid entity
//   cuboidEntity->addComponent(cuboidMesh);
//   cuboidEntity->addComponent(cuboidTransform);
//   cuboidEntity->addComponent(cuboidMaterial);

//   // Camera
//   Qt3DRender::QCamera *camera = view.camera();
//   camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f,
//   1000.0f); camera->setPosition(QVector3D(0, 0, 10));  // Initial camera
//   position camera->setViewCenter(QVector3D(0, 0, 0)); // Look at the origin

//   // Light
//   Qt3DCore::QEntity *lightEntity = new Qt3DCore::QEntity(rootEntity);
//   Qt3DRender::QPointLight *pointLight =
//       new Qt3DRender::QPointLight(lightEntity);
//   pointLight->setColor(Qt::white);
//   pointLight->setIntensity(10.0f);
//   lightEntity->addComponent(pointLight);

//   // Set the root entity of the scene
//   view.setRootEntity(rootEntity);

//   // Animation to rotate the camera around the cube
//   QPropertyAnimation *cameraAnimation =
//       new QPropertyAnimation(camera, "position");
//   cameraAnimation->setDuration(10000); // 10 seconds for a full rotation
//   cameraAnimation->setLoopCount(-1);   // Infinite loop

//   // Define the circular path for the camera
//   QVector3D center(0, 0, 0); // Center of the cube
//   float radius = 10.0f;      // Radius of the circular path
//   cameraAnimation->setKeyValueAt(
//       0, QVector3D(radius, 0, 0)); // Start at (radius, 0, 0)
//   cameraAnimation->setKeyValueAt(
//       0.25, QVector3D(10.0, radius, 10.0)); // Move to (0, radius, 0)
//   cameraAnimation->setKeyValueAt(
//       0.5, QVector3D(-radius, 0, 0)); // Move to (-radius, 0, 0)
//   cameraAnimation->setKeyValueAt(
//       0.75, QVector3D(0, -radius, 0)); // Move to (0, -radius, 0)
//   cameraAnimation->setKeyValueAt(
//       1, QVector3D(radius, 0, 0)); // Return to (radius, 0, 0)

//   // Start the animation
//   cameraAnimation->start();

//   // Show the window
//   view.show();

//   return app.exec();
// }
// //=========================END LINE
// DRAWING====================================
