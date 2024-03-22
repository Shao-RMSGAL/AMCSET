// main.cpp
#include <array>
#include <iostream>

extern "C" {
    void matrix_multiply(double a[6][6], double b[6][6], double c[6][6]);
    void initialize_identity(double a[6][6]);
}

int main() {
    double a[6][6], b[6][6], c[6][6];

    initialize_identity(a);
    initialize_identity(b);

    matrix_multiply(a, b, c);

    // Print c
    for (const auto& row : c) {
        for (const auto& elem : row) {
            std::cout << elem << ' ';
        }
        std::cout << '\n';
    }
}