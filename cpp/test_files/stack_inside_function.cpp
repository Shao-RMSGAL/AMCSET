#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int* foo() {
    int* ret = (int*)malloc(sizeof(int));
    *ret = 5;
    return ret;
}

int main() {
    int* j = foo();
    std::cout << *j << std::endl;
    free(j);
}