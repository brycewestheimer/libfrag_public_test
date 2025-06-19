#include "libfrag/libfrag.hpp"
#include "libfrag/libfrag_config.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
    std::cout << "libfrag Hello World Example" << std::endl;
    
    // Create an instance of MyClass
    libfrag::MyClass example(42);
    
    // Call the hello_world method
    example.hello_world();
    
    std::cout << "Example completed successfully!" << std::endl;
    
    return 0;
}