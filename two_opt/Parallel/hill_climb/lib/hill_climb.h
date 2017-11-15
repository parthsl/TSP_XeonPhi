#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits.h>
#include<stdbool.h>
#include<string.h>
#include<errno.h>
#include<float.h>
#include<omp.h>

/*
 * @Author : Parth Shah<parths1229@gmail.com>A
 *
 * Parallel VNN + 2opt(w/ Tiling and inline swapping) method combined method with OpenMP giving parallel execution of both in optimal way with co-ordinates stored for each city.
 * 2-opt works in tiles of data. Each thread has their own working set and tries to swap in that workspace without waiting for other threads.
 * 2-opt uses inline swapping technique : Swap all the edges giving the benefit.
 * Results : Speedup in VNN calculation. 2opt is more faster than w/o Tiling strategy + maxSwap technique but has increase in error rate. Increase in threads can give more faster result but increase in error from optimal distance too. Able to load and process 85k cities in few seconds.
 */

struct coords {
	double x,y;
};

struct coords* G;

typedef long int nd;

typedef volatile nd noopt;

#define euclidean_dist(i, j) \
        sqrt( (i.x-j.x)*(i.x-j.x) + (i.y-j.y)*(i.y-j.y))

void donot_optimize();

nd readGraph(char* filename,nd* cities);

double find_tour_length(struct coords* G , nd* tour, nd cities);

nd* VNN(struct coords* G, nd cities, nd start);

nd* VNNp(struct coords* G, nd cities, nd start);

void rev_arr(nd* min_circuit, nd s, nd e);

void print_tour(struct coords* G, nd* min_circuit, nd total_cities);

nd two_opt_inline_swap(struct coords* G, nd* min_circuit, nd cities);

nd two_opt_max_swap(struct coords* G, nd* min_circuit, nd cities);

void two_opt_random_swap(nd* min_circuit, nd cities, nd k);

void bruteforce(struct coords* G, nd* min_circuit, nd cities);

void chunk_bruteforce(struct coords* G, nd* min_circuit, nd cities, nd chunk_size);
