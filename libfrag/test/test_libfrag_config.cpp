#include <doctest.h>

#include "libfrag/libfrag.hpp"
#include "libfrag/libfrag_config.hpp"



TEST_SUITE_BEGIN("core");

TEST_CASE("check version"){

    #ifndef LIBFRAG_VERSION_MAJOR
        #error "LIBFRAG_VERSION_MAJOR is undefined"
    #endif
    

    #ifndef LIBFRAG_VERSION_MINOR
        #error "LIBFRAG_VERSION_MINOR is undefined"
    #endif


    #ifndef LIBFRAG_VERSION_PATCH
        #error "LIBFRAG_VERSION_PATCH is undefined"
    #endif

    CHECK_EQ(LIBFRAG_VERSION_MAJOR , 0);
    CHECK_EQ(LIBFRAG_VERSION_MINOR , 1);
    CHECK_EQ(LIBFRAG_VERSION_PATCH , 0);
}



TEST_SUITE_END(); // end of testsuite core
