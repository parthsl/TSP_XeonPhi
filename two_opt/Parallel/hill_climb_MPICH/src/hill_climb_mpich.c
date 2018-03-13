#include "hill_climb.h"

int main(int argc, char** argv) // argv1 = filename argv2 = threads
{
	int num_procs, num_local;
	nd total_cities = 0;
	nd* min_circuit = NULL;
	nd total_threads_available=1;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank (MPI_COMM_WORLD, &num_local);

	if(num_local==master)
		printf("Parallel ranks=%d\n",num_procs);

	srand(0);

	if(num_local == master && argc<2) {
		fprintf(stderr,"Provide input file\n");
		return 0;
	}
	if(num_local == master) {
		if( readGraph(argv[1],&total_cities)==EAGAIN) return EINVAL;
	}

	//Braodcast total cities in a GRAPH to every Node
	MPI_Bcast(&total_cities, 1, MPI_LONG_LONG, master, MPI_COMM_WORLD);

	if(num_local != master) {
		G = (double**)malloc(2*sizeof(double*));
		for(nd i=0; i<2; i++)
			G[i] = (double*)malloc(sizeof(double)*(total_cities+1));
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&G[0][0],total_cities+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&G[1][0],total_cities+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(argc>2)omp_set_num_threads(atoi(argv[2]));

	double stime = 0;
	if(num_local==master) stime = omp_get_wtime();


	#pragma omp parallel
	{
		total_threads_available = omp_get_num_threads();
	}
#ifdef DEBUG
	printf("%d: Threads spawn: %ld\n",num_local, total_threads_available);
#endif

	//---------------------------VNN---------------------------------------------------------------------//
	if(num_local==master) {
		omp_set_num_threads(8);
		min_circuit = VNN(G,total_cities,0);
		printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));
		omp_set_num_threads(total_threads_available);
	}
	else {
		min_circuit = (nd*)malloc(sizeof(nd)*(total_cities+1));
	}

	char mach_name[100];
	int mach_len;
	MPI_Get_processor_name(mach_name,&mach_len);
#ifdef DEBUG
	printf("%d:%s\n",num_local, mach_name);
#endif
	/*------------------------- Send min_circuit to everyone after VNN ---------------------------------*/
	MPI_Bcast(min_circuit,total_cities+1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

#ifdef DEBUG
	if(num_local==master)
		printf("Broadcasted Min_circuit id=%d\n",num_local);
#endif

	//---------------------------------------------------------------------------------------------------//
	MPI_Barrier(MPI_COMM_WORLD);

	if(num_local != master) {
		two_opt_max_swap(G, min_circuit,total_cities, num_procs, num_local);
	}
	else {
		nd counter = 0;
		double tour_length;
		counter += two_opt_max_swap(G, min_circuit,total_cities, num_procs, num_local);
		tour_length = find_tour_length(G,min_circuit,total_cities);
		printf("Hills climbed = %ld\n",counter);
		printf("Min distance = %lf\n",tour_length);
		printf("Time taken = %lf\n",omp_get_wtime()-stime);

#ifdef DEBUG
		print_tour(G,min_circuit,total_cities);
#endif
	}
	MPI_Finalize();
	free(min_circuit);
	free(G[0]);
	free(G[1]);
	return 0;
}
