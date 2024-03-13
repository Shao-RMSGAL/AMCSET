class B{
    private:
        int val;
    public:
        int getVal() { return val; };
};

class A{
    private:
        B* b;
    public:
        A(B* b) : b(b) {};
        static int foo(B* b) { return b->getVal(); };
};

int main() {
    B* beta = new B();
    A alpha(beta);
    alpha.foo(beta);

    return 0;
}