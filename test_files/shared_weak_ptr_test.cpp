#include <memory>
#include <iostream>

class Bombardment;

class Particle {
    
    private:
        std::weak_ptr<Bombardment> b_ref;
    
    public:
        Particle(std::shared_ptr<Bombardment> b_ref) : b_ref(b_ref) {};
        ~Particle() {
            std::cout << "Destroying Particle" << std::endl;
        };
};

class Bombardment : public std::enable_shared_from_this<Bombardment> {
    private:
        std::unique_ptr<Particle> a_ptr;
    
    public:
        Bombardment() {};
        ~Bombardment() {
            std::cout << "Destroying Bombardment" << std::endl;
        };

        void addParticle() {
            std::cout << "Adding particle" << std::endl;
             a_ptr = std::make_unique<Particle>(shared_from_this()); 
        };

};

int main() {
    std::shared_ptr<Bombardment> sharedB(new Bombardment());
    sharedB->addParticle();
    sharedB.reset();

    std::cout << "Test" << std::endl;
}