set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

add_executable(m_serial m_serial.c)
target_include_directories(m_serial PUBLIC ${INC} ${CMAKE_INSTALL_PREFIX})
target_link_libraries(m_serial m mandelbrot)

add_executable(m_static m_static.c)
target_include_directories(m_static PUBLIC ${INC} ${MPI_C_INCLUDE_PATH} ${CMAKE_INSTALL_PREFIX})
target_link_libraries(m_static ${MPI_C_LIBRARIES} m mandelbrot)