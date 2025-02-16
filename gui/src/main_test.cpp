#include <atomic>
#include <boost/asio.hpp>
#include <boost/asio/executor_work_guard.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <thread>
#include <variant>
#include <vector>

using DataVariant = std::variant<std::vector<int>, int, double>;
using namespace boost::asio;

struct ThreadTerminator {};
using QueueVariant = std::variant<DataVariant, ThreadTerminator>;

boost::lockfree::spsc_queue<QueueVariant *> queue(1024);
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

  void operator()(const ThreadTerminator &) const {
    std::cout << "Processing ThreadTerminator" << std::endl;
  }

  ~VariantVisitor() {
    std::cout << "Closing file" << std::endl;
    if (output_file.is_open()) {
      output_file.close();
    }
  }

private:
  void process_value(const std::vector<int> &vec) const {
    std::cout << "Writing vector to file" << std::endl;
    output_file << "Vector: ";
    for (int v : vec) {
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
    QueueVariant *raw_data = nullptr;
    if (queue.pop(raw_data)) {
      std::cout << "Popped item from queue, pending_items=" << pending_items
                << std::endl;
      std::unique_ptr<QueueVariant> data(raw_data);

      if (!std::holds_alternative<ThreadTerminator>(*data)) {
        std::cout << "Processing data item" << std::endl;
        std::visit(visitor, *data);
        pending_items--;
        std::cout << "Processed item, remaining pending_items=" << pending_items
                  << std::endl;
      }
    } else {
      std::this_thread::yield();
    }
  }
  std::cout << "Output thread terminating" << std::endl;
}

DataVariant bombardment_task() {
  std::cout << "Task starting" << std::endl;
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  std::cout << "Task finished" << std::endl;
  return std::vector<int>{1, 2, 3};
}

void simulation(io_context &io) {
  std::vector<std::thread> pool;
  const int num_threads = 4;
  const int tasks_per_thread = 10;

  // Use a work guard to keep the io_context active
  auto work = boost::asio::make_work_guard(io);

  std::cout << "Starting simulation threads" << std::endl;
  for (int i = 0; i < num_threads; ++i) {
    pool.emplace_back([&io, tasks_per_thread]() {
      for (int j = 0; j < tasks_per_thread; ++j) {
        std::cout << "Thread " << std::this_thread::get_id()
                  << " starting task " << j << std::endl;
        auto result = bombardment_task();
        post(io, [result]() {
          std::cout << "Post callback executing, incrementing pending_items"
                    << std::endl;
          pending_items++;
          auto *data = new QueueVariant(result);
          bool success = queue.push(data);
          if (!success) {
            std::cerr << "Failed to push to queue!" << std::endl;
            delete data;
            pending_items--;
          }
          std::cout << "Pushed to queue, pending_items=" << pending_items
                    << std::endl;
        });
      }
    });
  }

  // Wait for all pool threads to finish posting tasks
  std::cout << "Joining pool threads" << std::endl;
  for (auto &thread : pool) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  // Reset the work guard to allow the io_context to exit
  work.reset();

  std::cout << "Running io_context" << std::endl;
  io.run();

  std::cout << "All threads joined, sending terminator" << std::endl;
  should_terminate = true;
  auto *terminator = new QueueVariant(ThreadTerminator{});
  queue.push(terminator);
}

int main() {
  std::cout << "Program starting" << std::endl;
  io_context io;

  std::thread output_thread(output_thread_func);

  simulation(io);

  std::cout << "Waiting for output thread" << std::endl;
  if (output_thread.joinable()) {
    output_thread.join();
  }

  std::cout << "Program completed" << std::endl;
  return 0;
}
