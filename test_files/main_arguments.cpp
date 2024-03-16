#include <iostream>

int main(int argc, char* argv[]) {
    // Iterate through each command line argument
    for (int i = 0; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
    }

    return 0;
}