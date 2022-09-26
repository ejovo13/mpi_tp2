#include <stdio.h>
#include <stdint.h>
#include "ejovo.h"

#include <complex.h>
// I want to actually go ahead and use this file to draw the mandelbrot set.

// I think that I can easily get this work in within like 30 minutes
// I'll need complex numbers and access to my matrix functions

#define MAX_ITERATIONS 100
#define MAX_ABS_VALUE 20


// Simple iterate function
// We always start the iteration at z = 0. 
static inline double complex fc_z(double complex c, double complex z) {
    return z * z + c;
}

// Return the number of iterations for which c stays bouned i.e. cabs(fc_z(z) < MAX_ABS_VALUE)
int nb_iter(double complex c) {

    // We will always start iterating at z = 0
    double complex z = c;
    int nit = 0;

    while (cabs(z) < MAX_ABS_VALUE && nit < MAX_ITERATIONS) {
        z = fc_z(c, z); 
        nit ++;
    }

    return nit;
}

// Here we want to define a function that actually computes the grids coordinates
// given the width of the grid in pixels, the height of the grid in pixels, and the top left
// and bottom left points in C
Matrix_c *create_grid(double complex top_left, double complex bottom_right, size_t w_grid, size_t h_grid) {

    double cwidth  = creal(bottom_right) - creal(top_left);
    double cheight = cimag(top_left) - cimag(bottom_right);

    double dw =  cwidth / (w_grid - 1);
    double dh = cheight / (h_grid - 1);

    double x0 = creal(top_left), y0 = cimag(top_left);

    Matrix_c *grid = Matrix_new_c(h_grid, w_grid);

    for (size_t i = 0; i < grid->nrows; i++) {
        for (size_t j = 0; j < grid->ncols; j++) {
            matset_c( grid, i, j, CMPLX(x0 + j * dw, 
                                        y0 - i * dh) ); // since the imaginary part is decreasing, we subtract i
        }
    }

    return grid;
}

typedef int (*fn_toi_c) (double complex);

static inline Matrix_i *maptoi_c(const Matrix_c *A, fn_toi_c fn) {

    // allocate the new matrix
    Matrix_i *out = Matrix_new_i(A->nrows, A->ncols);
    const size_t n = Matrix_size_c(A);

    for (size_t i = 0; i < n; i++) {
        out->data[i] = fn(A->data[i]);
    }

    return out;
}


// Then we will need a different function to determine if a given c is in the mandelbrot set



int main() {

    // Let's define some numbers

    const double complex TOP_L = CMPLX(-2, 1.3);
    const double complex BOT_R = CMPLX(1, -1.3);

    const size_t PIXEL_H = 5000;
    const size_t PIXEL_W = PIXEL_H * 1.3;

    Matrix_c *grid = create_grid(TOP_L, BOT_R, PIXEL_W, PIXEL_H);
    // Map the nb_iter functions to the coordinate
    Matrix_i *iter_counts = maptoi_c(grid, nb_iter);

    write_ppm_grayscale(iter_counts, "mandelbrot5k.ppm", 255.0 / MAX_ITERATIONS);

    return 0;
}