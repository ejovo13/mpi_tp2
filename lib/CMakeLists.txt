add_library(mandelbrot mandelbrot_static.c mandelbrot.c)

target_include_directories(mandelbrot PUBLIC ${INC})
target_link_libraries(mandelbrot ${EJOVO_LIB})
