#ifndef MANDELBROT_H
#define MANDELBROT_H

#include <stdio.h>
#include <stdint.h>
#include "ejovo.h"

#include <complex.h>

#ifndef MAX_ITERATIONS
    #define MAX_ITERATIONS 200
#endif
#ifndef MAX_ABS_VALUE
    #define MAX_ABS_VALUE 2
#endif

typedef double complex c64;

// Simple iterate function
// We always start the iteration at z = 0. 
static inline c64 fc_z(c64 c, c64 z) {
    return z * z + c;
}

static inline int nb_iter_full(c64 c, double max_abs_value, int max_iterations) {

    // We will always start iterating at z = 0
    c64 z = c;
    int nit = 0;

    while (cabs(z) < max_abs_value && nit < max_iterations) {
        z = fc_z(c, z); 
        nit ++;
    }

    return nit;
}

// Return the number of iterations for which c stays bouned i.e. cabs(fc_z(z) < MAX_ABS_VALUE)
static inline int nb_iter(c64 c) {
    return nb_iter_full(c, MAX_ABS_VALUE, MAX_ITERATIONS);
}

// convert a linear idex to the row index (i) using a row_major indexing scheme
static inline int linear_index_to_i(int index, int grid_width) {
    return index / grid_width; //  
}

static inline int linear_index_to_j(int index, int grid_width) {
    return index % grid_width;
}

// Convert i, j coordinates to the corresponding complex number
static inline c64 grid_coords_to_complex(int i, int j, int grid_width, int grid_height, c64 top_left, c64 bottom_right) {

    double cwidth  = creal(bottom_right) - creal(top_left);
    double cheight = cimag(top_left) - cimag(bottom_right);

    double dw = cwidth / (grid_width - 1);
    double dh = cheight / (grid_height - 1);

    double x0 = creal(top_left), y0 = cimag(top_left);

    return CMPLX(x0 + j * dw, y0 - i * dh); // since the imaginary part is decreasing, we subtract i
}

// Add a function that converts a linear index (row-major) into a complex number
static inline c64 linear_index_to_complex(int index, int grid_width, int grid_height, c64 top_left, c64 bottom_right) {

    int i = linear_index_to_i(index, grid_width);
    int j = linear_index_to_j(index, grid_width);

    return grid_coords_to_complex(i, j, grid_width, grid_height, top_left, bottom_right);

}

static inline int count_iterations(int index, 
                                   int grid_width,
                                   int grid_height,
                                   c64 top_left, 
                                   c64 bottom_right,
                                   double max_abs_value,
                                   int max_iterations) 
{
    // Convert the linear index to a complex number
    c64 c = linear_index_to_complex(index, grid_width, grid_height, top_left, bottom_right);
    return nb_iter_full(c, max_abs_value, max_iterations);
}


// Here we want to define a function that actually computes the grids coordinates
// given the width of the grid in pixels, the height of the grid in pixels, and the top left
// and bottom left points in C
Matrix_c *create_grid(c64 top_left, c64 bottom_right, size_t w_grid, size_t h_grid);

typedef int (*fn_toi_c) (c64);

static inline Matrix_i *maptoi_c(const Matrix_c *A, fn_toi_c fn) {

    // allocate the new matrix
    Matrix_i *out = Matrix_new_i(A->nrows, A->ncols);
    const size_t n = Matrix_size_c(A);

    for (size_t i = 0; i < n; i++) {
        out->data[i] = fn(A->data[i]);
    }

    return out;
}

#endif // MANDELBROT_H