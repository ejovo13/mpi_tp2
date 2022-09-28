/**========================================================================
 * ?                          m_dynamic.c
 * @brief   : Test driver for dynamic load balancing
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-09-28
 *========================================================================**/

#include "mandelbrot_dynamic.h"
#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv){

    /**========================================================================
     *!                           MPI Initialization
     *========================================================================**/
    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    /**========================================================================
     *!                           Mandelbrot parameters
     *========================================================================**/
    const size_t GRID_H = 500;
    const size_t GRID_W = GRID_H * 1.3;
    const c64 TOP_L = CMPLX(-2, 1.3);
    const c64 BOT_R = CMPLX(0.5, -1.3);

    if (this_rank == 0) {
        manager(world_size, GRID_H, GRID_W);
    } else {
        worker(world_size, this_rank, GRID_H, GRID_W, TOP_L, BOT_R);
    }

    MPI_Finalize();
    return 0;
}