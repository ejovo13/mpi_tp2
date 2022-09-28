# Let's try using my ejovo library
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(Test_ejovo)

# set(CMAKE_PREFIX_PATH /usr/local/include/ejovo;/usr/local/include/ejovo/matrix)

set(INC ${CMAKE_SOURCE_DIR}/inc)
set(BIN ${CMAKE_BINARY_DIR})

# find_library(EJOVO_LIB ejovo)
find_package(MPI REQUIRED)

add_subdirectory(lib)
add_subdirectory(src)