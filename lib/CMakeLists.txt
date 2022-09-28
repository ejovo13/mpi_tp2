if (NOT EJOVO_LIB)
    message("ejovo not found, including libejovo subdirectory")
    add_subdirectory(libejovo)
endif()

add_library(mandelbrot mandelbrot_static.c mandelbrot.c)

target_include_directories(mandelbrot PUBLIC ${INC})

if (NOT EJOVO_LIB)
    target_link_libraries(mandelbrot ejovo)
else()
    target_link_libraries(mandelbrot ${EJOVO_LIB})
endif()
