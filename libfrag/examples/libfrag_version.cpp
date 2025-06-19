#include <iostream>

#include "libfrag/libfrag.hpp"
#include "libfrag/libfrag_config.hpp"

int main(int argc, char *argv[]){
    std::cout<<"LIBFRAG_VERSION_MAJOR "<<LIBFRAG_VERSION_MAJOR<<"\n";
    std::cout<<"LIBFRAG_VERSION_MINOR "<<LIBFRAG_VERSION_MINOR<<"\n";
    std::cout<<"LIBFRAG_VERSION_PATCH "<<LIBFRAG_VERSION_PATCH<<"\n";
}