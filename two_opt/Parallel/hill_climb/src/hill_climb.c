#include "hill_climb.h"

int main(int argc, char** argv) // argv1 = filename argv2 = threads
{
	nd total_cities, nn_routes=1;
	srand(0);
	if(argc<2) {
		fprintf(stderr,"Provide input file\n");
		return 0;
	}
	if(readGraph(argv[1],&total_cities)==EAGAIN)return EINVAL;
	if(argc>2)omp_set_num_threads(atoi(argv[2]));
	if(argc>3)nn_routes = atoi(argv[3])>0?atoi(argv[3]):1;

	double stime = omp_get_wtime();

	nd* min_circuit;
	nd total_threads_available=1;
	#pragma omp parallel
	{
		total_threads_available = omp_get_num_threads();
	}

	printf("Thread spawned: %ld\n",total_threads_available);

	double tour_length = DBL_MAX;

	#pragma omp parallel for schedule(static)
	for(int i=0; i<(nn_routes*omp_get_num_threads()); i++) {
		nd* nn_tour  = VNN(G,total_cities,i%total_cities);
		double nn_tour_cost = find_tour_length(G,nn_tour,total_cities);
		if(nn_tour_cost<tour_length)
			#pragma omp critical
		{
			if(nn_tour_cost<tour_length) {
				tour_length = nn_tour_cost;
				min_circuit = nn_tour;
			}
		}
		else free(nn_tour);
	}

	printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));

//---------------------------------------------------------------------------------------------------//

	nd counter = 0;

	counter += two_opt_max_swap(G, min_circuit,total_cities);
	tour_length = find_tour_length(G,min_circuit,total_cities);


	printf("Hills climbed = %ld\n",counter);
	printf("Min distance = %lf\n",tour_length);
	printf("Time taken = %lf\n",omp_get_wtime()-stime);
#ifdef DEBUG
	print_tour(G,min_circuit,total_cities);
#endif
	free(min_circuit);
	free(G);
	return 0;
}
