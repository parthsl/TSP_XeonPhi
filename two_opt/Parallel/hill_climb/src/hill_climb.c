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
	double tour_length;
	min_circuit = VNNp(G,total_cities,0);
	printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));
	
//---------------------------------------------------------------------------------------------------//

	double last_tour_length = DBL_MAX;
	nd total_threads_available=1;
	#pragma omp parallel 
	{
		total_threads_available = omp_get_num_threads();
	}
	nd counter = 0;

	while(total_threads_available>=1) {
		while(1) {
			counter += two_opt_inline_swap(G, min_circuit,total_cities);
			tour_length = find_tour_length(G,min_circuit,total_cities);
			if(last_tour_length<=tour_length)break;
			last_tour_length = tour_length;
		}
		printf("Distance optimized to = %lf\n",tour_length);
		total_threads_available/=2;
		omp_set_num_threads(total_threads_available);
	}

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
