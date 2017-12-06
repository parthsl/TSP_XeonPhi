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

	nd total_threads_available=1;
	#pragma omp parallel 
	{
		#pragma omp single nowait
		total_threads_available = omp_get_num_threads();
	}
	printf("Threads Spawn : %ld\n",total_threads_available);
	nd counter = 0;

	nd* optimal_tour = (nd*)malloc(sizeof(nd)*total_cities);
	double last_tour_length = DBL_MAX;
	nd k = 10;
	while(k<total_cities){
		while(1) {
			counter += two_opt_inline_swap(G, min_circuit,total_cities);
			tour_length = find_tour_length(G,min_circuit,total_cities);

			if(last_tour_length > tour_length){
				last_tour_length = tour_length;
				memcpy(optimal_tour,min_circuit,sizeof(nd)*total_cities);
			}
			else break;
		}
		printf("last_tour_length=%lf\n",last_tour_length);
		
		two_opt_random_swap(min_circuit, total_cities, k);

		while(1) {
                        counter += two_opt_inline_swap(G, min_circuit,total_cities);
                        tour_length = find_tour_length(G,min_circuit,total_cities);

                        if(last_tour_length > tour_length){
                                last_tour_length = tour_length;
                        }
                        else break;
                }

		//random swap giving any result then keep else revert back
		if(find_tour_length(G,min_circuit,total_cities)<find_tour_length(G,optimal_tour,total_cities)){
			memcpy(optimal_tour,min_circuit,sizeof(nd)*total_cities);

		}
		else{
			memcpy(min_circuit,optimal_tour,sizeof(nd)*total_cities);
		}
		k *= 2;
	}

	printf("Hills climbed = %ld\n",counter);
	printf("Min distance = %lf\n",find_tour_length(G,min_circuit,total_cities));
	printf("Time taken = %lf\n",omp_get_wtime()-stime);
#ifdef DEBUG
	print_tour(G,min_circuit,total_cities);
#endif
	free(min_circuit);
	free(G);
	return 0;
}
