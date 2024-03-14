//=======================MULTIPLE HANLDERS============================
// #include <boost/signals2.hpp>
// #include <iostream>

// void print_args(float x, float y)
// {
//   std::cout << "The arguments are " << x << " and " << y << std::endl;
// }

// void print_sum(float x, float y)
// {
//   std::cout << "The sum is " << x + y << std::endl;
// }

// void print_product(float x, float y)
// {
//   std::cout << "The product is " << x * y << std::endl;
// }

// void print_difference(float x, float y)
// {
//   std::cout << "The difference is " << x - y << std::endl;
// }

// void print_quotient(float x, float y)
// {
//   std::cout << "The quotient is " << x / y << std::endl;
// }

// int main() {
//   boost::signals2::signal<void (float, float)> sig;

//   sig.connect(&print_args);
//   sig.connect(&print_sum);
//   sig.connect(&print_product);
//   sig.connect(&print_difference);
//   sig.connect(&print_quotient);

//   sig(5., 3.);
// }

//=================BEGINNER EXAMPLE=================================

// #include <boost/signals2.hpp>
// #include <iostream>

// struct HelloWorld
// {
//   void operator()() const
//   {
//     std::cout << "Hello, World!" << std::endl;
//   }
// };

// int main() {

//     boost::signals2::signal<void ()> sig;
//     HelloWorld hello;
//     sig.connect(hello);

//     while(true) {
        
//     }
  
// }

#include <boost/signals2.hpp>
#include <iostream>
#include <boost/asio.hpp>

void clearLine() {
    constexpr int width = 80;
    std::cout << "\r";
    for (size_t i = 0; i < width; ++i) {
        std::cout << " ";
    }
    std::cout << "\r";
    std::cout.flush();
}

void handle_signal(const boost::system::error_code& /*ec*/, int /*signal_number*/)
{
    clearLine();
    std::cout << "Ctrl+C received. Exiting..." << std::endl;
    exit(0);
}

void mainloop(int argc, char* argv[]) {
    std::cout << argc << " argument" << ((argc==1) ? " " : "s ") << "came through. Here they are: " << std::endl;
    for(int i = 0; i < argc; i++) {
        std::cout << "[" << i << "]\t" << argv[i] << std::endl;
    }

    // // Read from standard input and echo to standard output
    std::string line;
    getline(std::cin, line);
    std::cout << "You typed: " << line << std::endl;
    std::cout << "Hi mom!" << std::endl;

}

int main(int argc, char* argv[]) {
    // Create boost::asio io_context
    boost::asio::io_context io_context;

    // Create a signal set
    boost::asio::signal_set signals(io_context, SIGINT);

    // Asynchronously wait for the signal to occur
    signals.async_wait(handle_signal);  

    // Submit a function to the io_context.
    boost::asio::post(io_context, [argc, argv]() {
        mainloop(argc, argv);
    });

    // Start the io_context
    io_context.run_one();

    return 0;
}
