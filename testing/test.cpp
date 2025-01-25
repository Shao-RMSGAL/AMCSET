#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

struct quaternion {
  double r;
  double x;
  double y;
  double z;
};

// Returns the product of two quaternions and returns the result
quaternion operator*(quaternion a, const quaternion &b) {
  auto r_1 = a.r;
  auto r_2 = b.r;
  auto x_1 = a.x;
  auto y_1 = a.y;
  auto z_1 = a.z;
  auto x_2 = b.x;
  auto y_2 = b.y;
  auto z_2 = b.z;

  auto r_3 = r_1 * r_2 - x_1 * x_2 - y_1 * y_2 - z_1 * z_2;
  auto x_3 = r_1 * x_2 + r_2 * x_1 + y_1 * z_2 - z_1 * y_2;
  auto y_3 = r_1 * y_2 + r_2 * y_1 + z_1 * x_2 - x_1 * z_2;
  auto z_3 = r_1 * z_2 + r_2 * z_1 + x_1 * y_2 - y_1 * x_2;

  return quaternion{r_3, x_3, y_3, z_3};
}

// Returns the sum of two quaternions
quaternion operator+(quaternion a, const quaternion &b) {
  return quaternion{a.r + b.r, a.x + b.x, a.y + b.y, +a.z + b.z};
}

// Returns the scalar multiplication of a quaternion
quaternion operator*(double a, const quaternion &b) {
  return quaternion{a * b.r, a * b.x, a * b.y, a * b.z};
}

// Returns the scalar division of a quaternion
quaternion operator/(const quaternion &b, double a) {
  return quaternion{b.r / a, b.x / a, b.y / a, b.z / a};
}

double quaternion_magnitude(const quaternion &a) {
  return std::sqrt(a.r * a.r + a.x * a.x + a.y * a.y + a.z * a.z);
}

// Returns the conjugtae of the quaternion
quaternion quaternion_conjugate(const quaternion &a) {
  return quaternion{a.r, -a.x, -a.y, -a.z};
}

// Returns a normalized quaternion
quaternion quaternion_normalize(const quaternion &a) {
  if (iszero(quaternion_magnitude(a))) {
    throw std::logic_error("Attempt to normalize 0-length quaternion");
  }
  auto result = a / quaternion_magnitude(a);
  return result;
}

quaternion operator!(const quaternion &a) {
  return quaternion{a.r, -a.x, -a.y, -a.z};
}

// Prints a quaternion
std::string print_quaternion(const quaternion &a) {
  std::stringstream stream;
  stream << a.r << "," << a.x << "," << a.y << "," << a.z;
  return stream.str();
}

size_t sysrandom(void *dst, size_t dstlen) {
  char *buffer = reinterpret_cast<char *>(dst);
  std::ifstream stream("/dev/urandom",
                       std::ios_base::binary | std::ios_base::in);
  stream.read(buffer, dstlen);
  return dstlen;
}

int main() {
  std::cout << "Starting quaternion experiment..." << std::endl;
  std::uint_least32_t seed;
  sysrandom(&seed, sizeof(seed));
  std::mt19937 gen(seed);
  std::uniform_real_distribution real_distribution;
  std::ofstream file_stream("output.csv");
  file_stream << "scalar,x,y,z\n";

  for (int particle = 0; particle < 100; particle++) {
    quaternion particle_position{0.0, 0.0, 0.0, 1.0};
    quaternion particle_velocity =
        quaternion_normalize(quaternion{0.0, 0.0, 0.0, 1.0});

    double speed = 1.0;

    for (int i = 0; i < 1000; i++) {
      // Generate random deflection angles
      double alpha_change =
          2 * M_PI * (0.5 - real_distribution(gen)); // Azimuthal
      // std::cout << "Alpha tilt: " << alpha_change * 180.0 / (2 * M_PI)
      //           << " degrees" << std::endl;

      double beta_change = 0.05 * real_distribution(gen); // Altitude
      // std::cout << "Beta tilt: " << beta_change * 180.0 / (2 * M_PI)
      //           << " degrees\n";

      // Quaternion for azimuthal rotation around the local z-axis
      auto q_alpha =
          quaternion{std::cos(alpha_change / 2.0),
                     particle_velocity.x * std::sin(alpha_change / 2.0),
                     particle_velocity.y * std::sin(alpha_change / 2.0),
                     particle_velocity.z * std::sin(alpha_change / 2.0)};

      // std::cout << "Alpha quaternion: " << print_quaternion(q_alpha)
      //           << " degrees\n";

      quaternion rotation_axis;
      // Calculate axis perpendicular to the current velocity for altitude
      // rotation
      if (particle_velocity.z == 1.0) {
        rotation_axis = quaternion{0.0, 1.0, 0.0, 0.0};

      } else {
        rotation_axis = quaternion_normalize(
            quaternion{0.0, particle_velocity.y, -particle_velocity.x, 0.0});
      }

      // std::cout << "Perpendicular axis: " << print_quaternion(rotation_axis)
      //           << " degrees\n";

      // Quaternion for altitude rotation around the perpendicular axis
      auto q_beta = quaternion{std::cos(beta_change / 2.0),
                               rotation_axis.x * std::sin(beta_change / 2.0),
                               rotation_axis.y * std::sin(beta_change / 2.0),
                               rotation_axis.z * std::sin(beta_change / 2.0)};

      // std::cout << "Beta quaternion: " << print_quaternion(q_beta)
      //           << " degrees\n";

      // Apply altitude rotation
      particle_velocity =
          q_alpha *
          ((q_beta * particle_velocity) * quaternion_conjugate(q_beta)) *
          quaternion_conjugate(q_alpha);

      particle_velocity.r = 0;
      particle_velocity = quaternion_normalize(particle_velocity);

      // Update particle position
      particle_position = particle_position + speed * particle_velocity;

      // Write results to the file
      file_stream << print_quaternion(particle_position) << "\n";
    }
  }

  std::cout << "...Calculation complete." << std::endl;

  return 0;
}
