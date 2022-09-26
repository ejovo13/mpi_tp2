/**========================================================================
 * ?                          mandelbrot.c
 * @brief   : Simple functions needed to draw the mandelbrot set
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-09-27
 *========================================================================**/
#include "mandelbrot.h"

Matrix_c *create_grid(c64 top_left, c64 bottom_right, size_t w_grid, size_t h_grid) {

    Matrix_c *grid = Matrix_new_c(h_grid, w_grid);

    for (size_t i = 0; i < grid->nrows; i++) {
        for (size_t j = 0; j < grid->ncols; j++) {
            matset_c( grid, i, j, grid_coords_to_complex(i, j, w_grid, h_grid, top_left, bottom_right) ); // since the imaginary part is decreasing, we subtract i
        }        
    }

    return grid;
}