#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits.h>
#include<stdbool.h>
#include<string.h>
#include<errno.h>
#include<float.h>
#include<omp.h>
#include<mpi.h>
#include<unistd.h>

#if defined(__ibmxl__) || defined(__powerpc__)
#include<massv.h>
#include<mass_simd.h>
#elif defined(__INTEL_COMPILER)
#include<mkl.h>
#endif

#define master 0

/*
 * @Author : Parth Shah<parths1229@gmail.com>A
 *
 * Parallel VNN + 2opt(w/ Tiling and inline swapping) method combined method with OpenMP giving parallel execution of both in optimal way with co-ordinates stored for each city.
 * 2-opt works in tiles of data. Each thread has their own working set and tries to swap in that workspace without waiting for other threads.
 * 2-opt uses inline swapping technique : Swap all the edges giving the benefit.
 * Results : Speedup in VNN calculation. 2opt is more faster than w/o Tiling strategy + maxSwap technique but has increase in error rate. Increase in threads can give more faster result but increase in error from optimal distance too. Able to load and process 85k cities in few seconds.
 */

double ** G;

typedef long int nd;

typedef volatile nd noopt;

#define squared_dist(G, i, j) \
        ( (G[0][i]-G[0][j])*(G[0][i]-G[0][j]) + (G[1][i]-G[1][j])*(G[1][i]-G[1][j]) )
#define euclidean_dist(G, i, j) \
        sqrt( (G[0][i]-G[0][j])*(G[0][i]-G[0][j]) + (G[1][i]-G[1][j])*(G[1][i]-G[1][j]) )

void donot_optimize();

nd readGraph(char* filename,nd* cities);

double find_tour_length(double** G , nd* tour, nd cities);

nd* VNN(double** G, nd cities, nd start);

nd* VNNp(double** G,nd cities, nd start);

void rev_arr(nd* min_circuit, nd s, nd e);

void print_tour(double** G, nd* min_circuit, nd total_cities);

nd two_opt_inline_swap(double** G, nd* min_circuit, nd cities);

nd two_opt_max_swap(double** G, nd* min_circuit, nd cities, int num_procs, int num_local);
