#include "hill_climb.h"

int main(int argc, char** argv) // argv1 = filename argv2 = threads
{
	nd total_cities;
	srand(0);
	if(argc<2) {
		fprintf(stderr,"Provide input file\n");
		return 0;
	}
	if(readGraph(argv[1],&total_cities)==EAGAIN)return EINVAL;
	if(argc>2)omp_set_num_threads(atoi(argv[2]));

	double stime = omp_get_wtime();

	nd* min_circuit;
	nd total_threads_available=1;
	printf("Thread spawned: %ld\n",total_threads_available);
  
	double tour_length;

	omp_set_num_threads(8);
	min_circuit = VNN(G,total_cities,0);
	omp_set_num_threads(total_threads_available);

	printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));
	
//---------------------------------------------------------------------------------------------------//

	nd counter = 0;
	
	counter += two_opt_max_swap_single(G, min_circuit,total_cities);
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
