# Let's try using my ejovo library


cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(Test_ejovo)

set(CMAKE_PREFIX_PATH /usr/local/include/ejovo;/usr/local/include/ejovo/matrix)

find_library(EJOVO_LIB ejovo REQUIRED)

message("Found lib: ${EJOVO_LIB}")

add_executable(test test.c)
target_link_libraries(test ${EJOVO_LIB} m)
target_include_directories(test PUBLIC ${CMAKE_PREFIX_PATH})