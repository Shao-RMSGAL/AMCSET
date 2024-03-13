#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

int main() {
    // Example string representing a file path
    std::string filePathString = "/path/to/file.txt";

    // Construct an std::filesystem::path object from the string
    fs::path filePath(filePathString);

    // Now you can use filePath as an std::filesystem::path object
    std::cout << "File name: " << filePath.filename() << std::endl;

    return 0;
}
