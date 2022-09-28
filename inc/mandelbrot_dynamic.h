/**========================================================================
 * ?                          mandelbrot_dynamic.h
 * @brief   : Dynamic load balancing to draw the mandelbrot set
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-09-28
 *========================================================================**/
#include <complex.h>
#include <mpi.h>

#include "mandelbrot.h"


#define MANAGER_WORK_ORDER 1
#define MANAGER_STOP_ORDER 2
// #define MANAGER_WORK_CONFIRMATION 3
// #define MANAGER_WORK_ORDER 2

// #define DEFAULT_GRID_H 5000
// #define DEFAULT_GRID_W 5000 * 1.3

void manager(int world_size, int const grid_h, int const grid_w) {

    MPI_Status status;
    int num_sent = 0;
    int sender, row;

    // Allocate space for the final image
    Matrix_i *img = Matrix_new_i(grid_h, grid_w);
    Matrix_i *img_row = Matrix_new_i(1, grid_w);
    
    // Loop through the rows, and send out the initial row number 
    for (int i = 1; i < world_size; i++) {
        MPI_Send(&i, 1, MPI_INT, i, MANAGER_WORK_ORDER, MPI_COMM_WORLD);
        num_sent++;
    }

    // Now receive the ALL the calculations!!
    for (int i = 0; i < grid_h; i++) {
        // Each iteration represents a SINGLE message exchange. Here we receive and later
        // in the loop we send another task to whichever process finished their calculation
        MPI_Recv(img_row->data, grid_w, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        sender = status.MPI_SOURCE;
        row = status.MPI_TAG;

        // Go ahead and fill the proper row of img
        for (int j = 0; j < grid_w; j++) {
            Matrix_set_i(img, row, j, img_row->data[j]);
        }

        // Now send a new task to the sender that finished
        if (num_sent < grid_h) {

            MPI_Send(&num_sent, 1, MPI_INT, sender, MANAGER_WORK_ORDER, MPI_COMM_WORLD); 
            num_sent++;

        } else { // no more work!!

            MPI_Send(MPI_BOTTOM, 0, MPI_INT, sender, MANAGER_STOP_ORDER, MPI_COMM_WORLD);
        }
    }

    write_ppm_grayscale(img, "m_dynamic.ppm", 255.0 / MAX_ITERATIONS);
}

void worker(int world_size, int this_rank, const int grid_h, const int grid_w, const c64 top_left, const c64 bottom_right) {

    // Each worker is going to receive the ROW that it needs to compute the number of iterations for
    int row;
    int rows_computed = 0;
    Matrix_i *img_row = Matrix_new_i(1, grid_w);
    MPI_Status status;

    // 
    if (this_rank <= grid_h) {

        // Receive an order
        MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        while (status.MPI_TAG == MANAGER_WORK_ORDER) {
            // Do some work

            // Compute the number of iterations for this row
            for (int j = 0; j < grid_w; j++) {

                c64 z = grid_coords_to_complex(row, j, grid_w, grid_h, top_left, bottom_right);
                img_row->data[j] = nb_iter(z);
            }

            // Send your results
            MPI_Send(img_row->data, grid_w, MPI_INT, 0, row, MPI_COMM_WORLD);

            rows_computed++;

            // And receive the next order
            MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        }
    }

    printf("p%d computed %d rows\n", this_rank, rows_computed);
}