/**========================================================================
 * ?                          mandelbrot_static.h
 * @brief   : Functions defining static load balancing used to compute
 *            a frame of the Mandelbrot set using MPI. 
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-09-27
 *========================================================================**/


#include "mandelbrot.h"

/**
 * @brief Compute the workload that `this_rank` will be responsible for
 * 
 * @param total_work Number of iterations that the entire program must 
 *                   divide amongst its processes
 * @param world_size total size of the MPI world
 * @param this_rank rank of the current MPI process
 * @return int The statically alloted workload to `this_rank`
 */
int compute_workload(int total_work, int world_size, int this_rank);

/**
 * @brief Compute an array that contains the workloads of each individual MPI process
 * 
 * @param total_work Total work (in number of iterations) that this program must execute
 * @param world_size size of the MPI world
 * @return Matrix_i* workloads of each process, in order 
 */
Matrix_i *compute_workload_array(int total_work, const int world_size);

/**
 * @brief Compute the displacements array needed by `MPI_Gatherv`
 * 
 * @param total_work Total work (in number of iterations) that this program must execute
 * @param world_size size of the MPI world
 * @return Matrix_i* Integer array of displacements. Be mindful to select the actual `data` 
 *         field of this matrix when passing the result to MPI_Gatherv
 */
Matrix_i *compute_displacements(int total_work, const int world_size);

/**
 * @brief Compute an array that simultaneously contains the start indices and end indices of all the processes.
 * 
 * The start index of a process p can be retrieved by accessing the returned matrix's pth element + 1; The end index
 * of a process p can be retrieved by accessing the returned matrix's (pth + 1) element.
 * 
 * @example 
 * ```
 * Matrix_i *startend = compute_startend_array(500, 7);
 * int start = startend->data[p] + 1;
 * int end   = startend->data[p + 1];
 * ``` 
 * 
 * @param total_work 
 * @param world_size 
 * @return Matrix_i* 
 */
Matrix_i *compute_startend_array(int total_work, int world_size);
