#include <memory>
#include <iostream>
#include <stdlib.h>

int main() {
    char* test_array = (char*)malloc(10);

    free(test_array);
}