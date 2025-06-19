
cmake_minimum_required(VERSION 3.15)  # Modern CMake version
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(libfrag LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)  # More modern than C++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project directories
set(PROJECT_INCLUDE_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CPP_ROOT_FOLDER_NAME   libfrag)
include_directories(${PROJECT_INCLUDE_DIR})   


# Versionning
# ===========
set(CPP_CONFIG_HEADER_NAME libfrag_config.hpp)
set(PROJECT_CONFIG_FILE ${PROJECT_INCLUDE_DIR}/${CPP_ROOT_FOLDER_NAME}/${CPP_CONFIG_HEADER_NAME})

file(STRINGS ${PROJECT_CONFIG_FILE} project_version_defines
     REGEX "#define LIBFRAG_VERSION_(MAJOR|MINOR|PATCH)")
message(STATUS ${PROJECT_CONFIG_FILE})

foreach(ver ${project_version_defines})
    if(ver MATCHES "#define LIBFRAG_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(${PROJECT_NAME}_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()

set(${PROJECT_NAME}_VERSION
    ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH})
message(STATUS "Building project v${${PROJECT_NAME}_VERSION}")

# Dependencies
# ============
find_package(xtl      REQUIRED)
find_package(xtensor  REQUIRED)

# # Build
# # =====

SET(INTERFACE_LIB_NAME libfrag)

# Collect headers and source files
file(GLOB_RECURSE ${PROJECT_NAME}_HEADERS  ${PROJECT_INCLUDE_DIR}/*.hpp)
file(GLOB_RECURSE ${PROJECT_NAME}_SOURCES  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Create library - now a regular library with implementation files
add_library(${INTERFACE_LIB_NAME} ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})

message(STATUS "INSTALL_INTERFACE: ${CMAKE_INSTALL_INCLUDEDIR}")
target_include_directories(${INTERFACE_LIB_NAME} PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Set C++ standard for the library
target_compile_features(${INTERFACE_LIB_NAME} PUBLIC cxx_std_17)

target_link_libraries(${INTERFACE_LIB_NAME} 
  PUBLIC xtensor)



# Benchmark
# ============
OPTION(BUILD_DOCS "build docs" ON)

OPTION(BUILD_BENCHMARK "${PROJECT_NAME} test suite" ON)
OPTION(DOWNLOAD_GOOGLE_BENCHMARK "build google benchmark from downloaded sources" ON)

if(DOWNLOAD_GOOGLE_BENCHMARK)
    set(BUILD_BENCHMARK ON)
endif()

if(BUILD_BENCHMARK)
    add_subdirectory(benchmark)
endif()

if(BUILD_DOCS)
    add_subdirectory(docs)
endif()

# Python
# ============
OPTION(BUILD_PYTHON "${PROJECT_NAME} python binding" ON)
if(BUILD_PYTHON)
    add_subdirectory(python)
endif()



# Tests
###########

OPTION(BUILD_TESTS "${PROJECT_NAME} test suite" ON)
# OPTION(DOWNLOAD_DOCTEST "build doctest from downloaded sources" ON)

# if(DOWNLOAD_DOCTEST)
#     set(BUILD_TESTS ON)
# endif()

if(BUILD_TESTS)
    add_subdirectory(test)
endif()

# Examples
# ============
OPTION(BUILD_EXAMPLES "${PROJECT_NAME} test suite" ON)


if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()



# Installation
# ============

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


set(${PROJECT_NAME}_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for ${PROJECT_NAME}Config.cmake")


message(STATUS "DEST ${${PROJECT_NAME}_CMAKECONFIG_INSTALL_DIR}" )


install(TARGETS ${INTERFACE_LIB_NAME}
        EXPORT ${INTERFACE_LIB_NAME}-targets)

install(EXPORT ${INTERFACE_LIB_NAME}-targets
        FILE ${INTERFACE_LIB_NAME}Targets.cmake
        #NAMESPACE ${PROJECT_NAME}:: 
        DESTINATION lib/cmake/${PROJECT_NAME}
         )


install(DIRECTORY  ${PROJECT_INCLUDE_DIR}/${CPP_ROOT_FOLDER_NAME}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})




configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${${PROJECT_NAME}_CMAKECONFIG_INSTALL_DIR})


# ${PROJECT_NAME} is header-only and does not depend on the architecture.
# Remove CMAKE_SIZEOF_VOID_P from ${PROJECT_NAME}ConfigVersion.cmake so that an ${PROJECT_NAME}Config.cmake
# generated for a 64 bit target can be used for 32 bit targets and vice versa.
set(_${PROJECT_NAME}_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                 VERSION ${${PROJECT_NAME}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
set(CMAKE_SIZEOF_VOID_P ${_${PROJECT_NAME}_CMAKE_SIZEOF_VOID_P})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${${PROJECT_NAME}_CMAKECONFIG_INSTALL_DIR})

install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${${PROJECT_NAME}_CMAKECONFIG_INSTALL_DIR})

configure_file(${PROJECT_NAME}.pc.in
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
                @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")



