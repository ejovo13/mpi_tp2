/**========================================================================
 * ?                          mandelbrot_static.c
 * @brief   : Parallel Mandelbrot set using static load-balancing
 * @details : This program utilizes MPI to create one frame of the mandelbrot 
 *            set. 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-09-26
 *========================================================================**/
#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include <mpi.h>

#include "mandelbrot_static.h"

int main(int argc, char** argv) {

    /**========================================================================
     *!                           Grid Parameters
     *========================================================================**/
    const size_t GRID_H = 5000;
    const size_t GRID_W = 1.3 * GRID_H;
    const c64 TOP_L = CMPLX(-2, 1.3);
    const c64 BOT_R = CMPLX(1, -1.3);

    const size_t TOTAL_WORK = GRID_H * GRID_W;

    /**========================================================================
     *!                           MPI setup
     *========================================================================**/
    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    Matrix_i *startend = compute_startend_array(TOTAL_WORK, world_size);

    // Really tough lesson to learn. If the matrices haven't been allocated and we are working with null pointers 
    // then we can't actually access their addresses. So, if we don't instantiate these objects
    // then we will end up with an error.
    Matrix_i *gathered_iterations = Matrix_new_i(0, 0);
    Matrix_i *recv_counts = Matrix_new_i(0, 0);
    Matrix_i *displacements = Matrix_new_i(0, 0); // ONLY TO BE ALLOCATED BY RANK 0 FOR GATHER OPERATION

    /**========================================================================
     *!                          Static Load-balancing
     *========================================================================**/
    // Allocate the appropriate amount of space for the matrices that are needed by MPI_Gatherv
    if (this_rank == 0) {

        Matrix_free_i(gathered_iterations);
        Matrix_free_i(recv_counts);
        Matrix_free_i(displacements);

        gathered_iterations = Matrix_new_i(GRID_H, GRID_W);
        displacements = compute_displacements(TOTAL_WORK, world_size);   
        recv_counts = compute_workload_array(TOTAL_WORK, world_size);    

    }

    MPI_Barrier(MPI_COMM_WORLD);

    int start = startend->data[this_rank] + 1;
    int end   = startend->data[this_rank + 1];
    int workload = (end - start) + 1;

    for (int i = 0; i < world_size; i++) {

        if (i == this_rank) {
            fprintf(stderr, "p%d processing [%d, %d]\t (%d elements)\n", this_rank, start, end, workload);
        };

        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (this_rank == 0) {
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /**========================================================================
     *!                           Computation of iterations
     *========================================================================**/

    // Now that the workloads have been appropriately distributed, we want to go ahead and allocate the matrix
    // that will store the iterations for THIS workload.
    // Since I have a "workload" matrix, I actually don't have to fuss around with a matrix that stores indices
    Matrix_i *nb_iterations = Vector_new_i(workload);
    Clock *clock = Clock_new();

    // Using the mandelbrot functions, go ahead and compute the indices that were computed
    Clock_tic(clock);
    for (int i = 0, mandelbrot_index = start; i < workload; i++, mandelbrot_index++) {
        
        nb_iterations->data[i] = count_iterations(mandelbrot_index, 
                                                  GRID_W,
                                                  GRID_H,
                                                  TOP_L,
                                                  BOT_R,
                                                  MAX_ABS_VALUE,
                                                  MAX_ITERATIONS);
    }
    Clock_toc(clock);

    MPI_Barrier(MPI_COMM_WORLD);

    /**========================================================================
     *!                 Start cleaning up; write to .ppm
     *========================================================================**/
    for (int i = 0; i < world_size; i++) {

        if (i == this_rank) {
            fprintf(stderr, "p%d finished in %lf s\n", this_rank, elapsed_time(clock));
        };

        MPI_Barrier(MPI_COMM_WORLD);

    }

    // Since nb_iterations, gathered_iterations, recv_counts, and displacements are 
    // `Matrix_i` structures (defined in libejovo), we pass their underlying data arrays
    // in accordance with MPI_Gatherv's signature
    MPI_Gatherv(nb_iterations->data,
                workload,
                MPI_INT,
                gathered_iterations->data,
                recv_counts->data,
                displacements->data,
                MPI_INT,
                0,
                MPI_COMM_WORLD
                );
            
    // Now that I've gathered the data, let's go ahead and actually output it.
    if (this_rank == 0) {
        write_ppm_grayscale(gathered_iterations, "m_static.ppm", 255.0 / MAX_ITERATIONS);
    }

    MPI_Finalize();
    return 0;
}