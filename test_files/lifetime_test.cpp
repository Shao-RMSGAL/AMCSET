#include <iostream>

class A{
    private:
        int memberA;
    public:
        A(const int& memberA) : memberA(memberA) {};

        void foo() {
            int a = 3;
            updateA(a);
        };

        void updateA(const int& a) {
            memberA = a;
        };

        void print() const {
            std::cout << memberA << std::endl;
        };
};

int main() {
    A instanceA(5);
    instanceA.foo();
    instanceA.print();
}