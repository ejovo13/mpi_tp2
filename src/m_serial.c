/**========================================================================
 * ?                          mandelbrot_serial.c
 * @brief   : Serial version of the Mandelbrot set
 * @details : Driving program for the functions that are declared in
 *            mandelbrot.h. This program will compute the number of 
 *            iterations that it takes for a complex number to diverge
 *            (or not diverge) when we reapply the complex valued function
 *            f_c(z) = z**2 + c. We discretize a grid based on the values
 *            of TOP_L, BOT_R, PIXEL_H, and PIXEL_W; then we output the matrix
 *            of iterations to a .ppm file that can be converted to png using
 *            a variety of tools 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-09-26
 *========================================================================**/
#include <complex.h>

#include "mandelbrot.h"

int main() {

    const size_t GRID_H = 5000;
    const size_t GRID_W = GRID_H * 1.3;
    const c64 TOP_L = CMPLX(-2, 1.3);
    const c64 BOT_R = CMPLX(1, -2);
    const char filename[] = "m_serial.ppm";

    Clock *clock = Clock_new();

    printf("Creating Mandelbrot grid with resolution %lu x %lu\n", GRID_H, GRID_W);


    fprintf(stderr, "Discretizing grid... ");
    Clock_tic(clock); 
    Matrix_c *grid = create_grid(TOP_L, BOT_R, GRID_W, GRID_H);
    Clock_toc(clock);

    fprintf(stderr, "Done. (%lf s)\n", elapsed_time(clock)); 

    fprintf(stderr, "Mapping complex values to number of iterations... ");
    Clock_tic(clock);
    Matrix_i *iter_counts = maptoi_c(grid, nb_iter); // map nb_iter to the elements of grid
    Clock_toc(clock);

    fprintf(stderr, "Done. (%lf s)\n", elapsed_time(clock));

    fprintf(stderr, "Writing to '%s'... ", filename);
    Clock_tic(clock);
    write_ppm_grayscale(iter_counts, filename, 255.0 / MAX_ITERATIONS);
    Clock_toc(clock);
    fprintf(stderr, "Done! (%lf)\n", elapsed_time(clock));

    return 0;
}