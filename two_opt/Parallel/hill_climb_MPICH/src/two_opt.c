#include "hill_climb.h"

nd two_opt_inline_swap(double** G, nd* min_circuit, nd cities) {
	nd counter_global = 0;
	noopt* lock;
#pragma omp parallel
	{
		nd counter = 0;
		nd id = omp_get_thread_num();
		nd total_threads;
		nd block_size;

		total_threads = omp_get_num_threads();
		if(total_threads==0)
			block_size = 1;
		else block_size = cities/total_threads;

#pragma omp single
		lock = (noopt*)calloc(total_threads,sizeof(noopt));

		nd iend;
		if(id==total_threads-1)
			iend = cities+1;
		else iend = block_size*(id+1);
#pragma omp barrier

		nd i = block_size*id;
		double max_change_local = 0;

		for(; i<iend-2; i++) {
			nd i_city = min_circuit[i];
			for(nd j=i+2; j<iend-1; j++) {
				nd j_city = min_circuit[j];
				nd j_next_city = min_circuit[j+1];
				nd i_next_city = min_circuit[i+1];
				double s_dist = euclidean_dist(G, i_city,j_city) + euclidean_dist(G, i_next_city, j_next_city);
				double f_dist = euclidean_dist(G, i_city, i_next_city) + euclidean_dist(G, j_city, j_next_city);
				if(f_dist>s_dist) {
					if(f_dist-s_dist > max_change_local) {
						for(nd z=0; z<=(j-i-1)/2; z++) {
							nd temp = min_circuit[i+1+z];
							min_circuit[i+1+z] = min_circuit[j-z];
							min_circuit[j-z] = temp;
						}
						counter++;
					}
				}
			}
		}

		lock[id]=block_size*id;

		if(id<total_threads-1) {

			nd j = block_size*(id+1);

			for(nd g=block_size*id+1; g<iend; g++) { // Runs through whole block

				while(lock[id+1]<j);
				nd j_next_city = min_circuit[j+1];
				for(i=g; i<iend; i++) // runs i in whole block but taking j+1 as const
				{
					nd i_city = min_circuit[i];
					nd i_next_city = min_circuit[i+1];
					nd j_city = min_circuit[j];

					double s_dist = euclidean_dist(G, i_city, j_city) + euclidean_dist(G,i_next_city, j_next_city);
					double f_dist = euclidean_dist(G, i_city, i_next_city) + euclidean_dist(G, j_city, j_next_city);
					if(f_dist>s_dist) {
						if(f_dist-s_dist > max_change_local) {
							for(nd z=0; z<=(j-i-1)/2; z++) {
								nd temp = min_circuit[i+1+z];
								min_circuit[i+1+z] = min_circuit[j-z];
								min_circuit[j-z] = temp;
							}
							counter++;
						}
					}
				}
				j++;
				lock[id]++;
			}

		}//End IF
		else {
			nd j = 0;

			for(nd g=block_size*id+1; g<block_size*(id+1); g++) { // Runs through whole block
				while(lock[0]<j);
				nd j_next_city = min_circuit[j+1];
				for(i=g; i<iend-2; i++) // runs i in whole block but taking j+1 as const
				{
					nd i_city = min_circuit[i];
					nd i_next_city = min_circuit[i+1];
					nd j_city = min_circuit[j];

					double s_dist = euclidean_dist(G, i_city, j_city) + euclidean_dist(G,i_next_city, j_next_city);
					double f_dist = euclidean_dist(G, i_city, i_next_city) + euclidean_dist(G, j_city, j_next_city);
					if(f_dist>s_dist) {
						for(nd z=0; z<=(cities+j-i-1)/2; z++) {
							nd temp = min_circuit[(i+1+z)%cities];
							min_circuit[(i+1+z)%cities] = min_circuit[(cities+j-z)%cities];
							min_circuit[(cities+j-z)%cities] = temp;
						}
						counter++;
					}
				}
				j++;
				lock[id]++;
			}

		}

#pragma omp critical
		counter_global += counter;
	}//pragma over

	return counter_global;
}


nd two_opt_max_swap(double** G, nd* min_circuit, nd cities, int num_procs, int num_local) {
	//for Mpich
	nd alpha2 = (1- (float)(num_local)/num_procs)*(cities-2)*(cities-3);
	nd sub_block_start = cities - ceil(( sqrt((alpha2-6)*4 + 25) +5) /2);
	if(num_local == 0)sub_block_start=0;
	alpha2 = (1- (float)(num_local+1)/num_procs)*(cities-2)*(cities-3);
	nd sub_block_end = cities - ceil(( sqrt((alpha2-6)*4 + 25) +5) /2);
	if(num_local == num_procs-1){
		sub_block_end = cities-1;
	}

#if defined(DEBUG) && defined(DEBUGL1)
	printf("Rank %d: Process Working-Set Block Start\\End=%ld\\%ld\n",num_local,sub_block_start,sub_block_end);
#endif

	//for Openmp
	double max_change = 0;
	nd ic=0,jc=0;
	nd counter = 0;
	double precal_distance[cities];
	enum loopFlags{ Proceed, Abort};
	static int loop_flag= Proceed;

#pragma omp parallel
	{
		nd id = omp_get_thread_num();
		nd total_threads = omp_get_num_threads();
		const long long int VS = 64;
		int VVS = VS;
		double avx_ed[VS];
		double avx_pre[VS];

		//Find iend
		nd out_of_scope_iterations = (cities-sub_block_end-2)*(cities-sub_block_end-3);
		if( (cities -sub_block_end -3) < 0) out_of_scope_iterations = 0;
		nd block_iterations = (cities-sub_block_start-2)*(cities-sub_block_start-3) - out_of_scope_iterations;
		nd alpha = out_of_scope_iterations + (1 - (float)(id+1)/total_threads)*block_iterations;
		nd iend = cities - ceil(( sqrt((alpha-6)*4 + 25) +5) /2);
		if(id==total_threads-1) iend = sub_block_end;

#if defined(__ibmxl__) || defined(__powerpc__)
                vector double x = {0,0},y={0,0};
#endif

		//build precalculated euclidean distance to neighbours
		#pragma omp for
                for(nd i=0; i<cities; i++) precal_distance[i] = euclidean_dist(G, min_circuit[i],min_circuit[i+1]);
#pragma omp barrier
		//Find istart
		alpha = out_of_scope_iterations +(1 - (float)(id)/total_threads)*block_iterations;
		alpha = cities - ceil(( sqrt((alpha-6)*4 + 25)+5) /2);
		if(id==0)alpha=sub_block_start;

#if defined(DEBUG) && defined(DEBUGL1)
		printf("ThreadID:%ld iteration_start\\end:%ld\\%ld\n",id,alpha,iend);
#endif

		//Optimised while loop with avx2
		while(loop_flag != Abort) {
			nd ic_local=0, jc_local=0;
			nd i = alpha;
			double max_change_local = 0;
			
			for(; i<iend; i++) { //goto iend as OutOfRange Checking is done by inner j loop
				nd i_city = min_circuit[i];
				nd i_next_city = min_circuit[i+1];
				nd j = i+2;

				for(VVS=VS; VVS>=2; VVS=VVS/2) {
                                        for(; j<(cities-1-VVS); j+=VVS) {
#if defined(__ibmxl__) || defined(__powerpc__)
                                                for(nd jj=0; jj<VVS; jj++) {
                                                        avx_ed[jj] = squared_dist(G, i_city, min_circuit[j+jj]);
                                                        avx_pre[jj] = squared_dist(G, i_next_city, min_circuit[j+jj+1]);
                                                }
                                                vsqrt(avx_ed,avx_ed,&VVS);
                                                vsqrt(avx_pre,avx_pre,&VVS);
                                                for(nd jj=0; jj<VVS; jj++) avx_ed[jj] = avx_ed[jj] + avx_pre[jj];
#else
                                                for(nd jj=0; jj<VVS; jj++) {
                                                        avx_ed[jj] = euclidean_dist(G, i_city, min_circuit[j+jj]) + euclidean_dist(G, i_next_city, min_circuit[j+jj+1]);
                                                }
#endif

                                                for(nd jj=0; jj<VVS; jj++) {
                                                        avx_pre[jj] = precal_distance[i] + precal_distance[j+jj];
                                                }
                                                for(nd jj=0; jj<VVS; jj++) {
                                                        avx_pre[jj] = avx_pre[jj] - avx_ed[jj];
                                                }
                                                for(nd jj=0; jj<VVS; jj++) {
                                                        if(avx_pre[jj] > max_change_local) {
                                                                max_change_local = avx_pre[jj];
                                                                ic_local = i;
                                                                jc_local = j+jj;
                                                        }
                                                }
                                        }
                                }

				//Peeled loop
				for(; j<cities-1; j++) {
					nd j_city = min_circuit[j];
					nd j_next_city = min_circuit[j+1];
#if defined(__ibmxl__) || defined(__powerpc__)
                                        x[0] = squared_dist(G, i_city, j_city);
                                        x[1] = squared_dist(G, i_next_city, j_next_city);
                                        y = sqrtd2(x);
                                        double s_dist = y[0] + y[1];
#else
                                        double s_dist = euclidean_dist(G, i_city, j_city) + euclidean_dist(G, i_next_city, j_next_city);
#endif
					double f_dist = precal_distance[i] + precal_distance[j];
					if(f_dist>s_dist) {
						if(f_dist-s_dist > max_change_local) {
							max_change_local = f_dist-s_dist;
							ic_local = i;
							jc_local=j;
						}
					}
				}
			}

#pragma omp critical
			{
				if(max_change_local > max_change) {
					max_change = max_change_local;
					ic = ic_local;
					jc = jc_local;
				}
			}

#pragma omp barrier
			if(id==0)
				MPI_Barrier(MPI_COMM_WORLD);
#pragma omp barrier

/*
 * Only ID=0, i.e. only main root OpenMP thread will go in this section
 * MPI Communication takes place in this block to find best swap among all nodes
 * Works with only one Global Communicator
 * Protocol:
 * 	1. MASTER Node collects all the max_swap from every node and stores in a array
 *	2. MASTER Node finds a node with max_swap
 *	3. MASTER Node broadcasts Max_Change value possible with max_swap to every node.
 *	4. MASTER Node broadcasts Node selected for max_swap
 *	5. Choosen Node Sends its ic and jc co-ordinates of the edges giving max_swap to MASTER Node.
 *	6. MASTER Node receives this co-ordinates and Broadcasts this to every node.
 *	7. MASTER Node decides to further proceed or not based on Max_Change value and broadcasts to every node.
 *	8. Every Node recieves loop_flag from MASTER Node indicating to preceed further or abort the execution. On Proceed, they swapps the choosen edges by reversing the min_circuit sub-array and pre-calculate array to maintain coherency.
 *	9. Every Node on Abort stops the execution and comes out of Computation loop.
 */
			if(id==0){
			MPI_Barrier(MPI_COMM_WORLD);
				double* max_change_procs = NULL;
				nd selected_coords[2];
				nd selected_nodes_swap = 0;
				if(num_local == master)
					max_change_procs = (double*)malloc(num_procs*sizeof(double));

				MPI_Gather(&max_change, 1, MPI_DOUBLE, max_change_procs, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);///Step 1. Master Node will save all distance computed by other node


				if(num_local == master){ //Step 2. Find node with maxswap
					for(nd i = 0; i<num_procs; i++){
						if(max_change < max_change_procs[i]){
							selected_nodes_swap = i;
							max_change = max_change_procs[i];
						}
					}

					free(max_change_procs);
				}

				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Bcast(&max_change, 1, MPI_DOUBLE, master, MPI_COMM_WORLD); //Step 3. Send Max Change Distance
				MPI_Bcast(&selected_nodes_swap, 1, MPI_LONG_LONG, master, MPI_COMM_WORLD); //Step 4. Send Selected node in MPI Cluster giving best swap

				selected_coords[0] = ic;
				selected_coords[1] = jc;
				if(num_local == selected_nodes_swap && num_local != master ){ // Step 5. Selected Node sends ic and jc from Selected node to MASTER node
#if defined(DEBUG) && defined(MPI_DEBUG)
					printf("%d sends ic=%ld jd=%ld to master\n",num_local,selected_coords[0],selected_coords[1]);
#endif
					MPI_Send(selected_coords,2,MPI_LONG_LONG, master, 0, MPI_COMM_WORLD);
				}
				else if(num_local != selected_nodes_swap && num_local==master){ //Step 6. Master Node receives the selected ic and jc cities
					MPI_Recv(selected_coords, 2, MPI_LONG_LONG, selected_nodes_swap,0,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if defined(DEBUG) && defined(MPI_DEBUG)
					printf("Master recieves ic=%ld jc=%ld\n",selected_coords[0],selected_coords[1]);
#endif
				}

#if defined(DEBUG) && defined(MPI_DEBUG)
				if(num_local==master)
					printf("Broadcasting selected ic and jc to all other nodes\n");
#endif

				MPI_Bcast(selected_coords, 2, MPI_LONG_LONG, master, MPI_COMM_WORLD); //Step 6. Broadcast Selected ic and jc to each nodes
				ic = selected_coords[0];
				jc = selected_coords[1];
				if(num_local == master && max_change > 0)//Step 7. Master Node decides to further Proceed or not by setting loop_flag variable
					loop_flag = Proceed;
				else loop_flag = Abort;

#if defined(DEBUG) && defined(MPI_DEBUG)
				if(num_local==master)
					printf("Broadcasting loop flag=%d\n",loop_flag);
#endif

				MPI_Bcast(&loop_flag, 1, MPI_INT, master, MPI_COMM_WORLD); //Step 7. Broadcast whether to proceed or stop

				/*Step 8. Every Node swaps edges if MASTER told to Proceed*/

				if(loop_flag == Proceed){
					static nd j;
					j = (jc-ic-1)/2;
					for(i=0; i<=j; i++) {
						nd temp = min_circuit[ic+1+i];
						min_circuit[ic+1+i] = min_circuit[jc-i];
						min_circuit[jc-i] = temp;
					}
				 for(i=0; i<=j; i++) {
	                                        nd temp = precal_distance[ic+i];
        	                                precal_distance[ic+i] = precal_distance[jc-i];
                	                        precal_distance[jc-i] = temp;
                        	        }
#if defined(__ibmxl__) || defined(__powerpc__)
                                	x[0] = squared_dist(G, min_circuit[ic], min_circuit[ic+1]);
	                                x[1] = squared_dist(G, min_circuit[jc], min_circuit[jc+1]);
        	                        y = sqrtd2(x);
                	                precal_distance[ic] = y[0];
                        	        precal_distance[jc] = y[1];
#else
                                	precal_distance[ic] = euclidean_dist(G, min_circuit[ic], min_circuit[ic+1]);
	                                precal_distance[jc] = euclidean_dist(G, min_circuit[jc], min_circuit[jc+1]);
#endif

					max_change = 0;
					ic=0;jc=0;
					counter++;
				}
				else{ //Step 9. MASTER told to abort the execution and terminate all threads as no max_swap benefit obtained now.
					loop_flag = Abort;
				}
			}//if over for id==0
#pragma omp barrier
			if(id==0)
				MPI_Barrier(MPI_COMM_WORLD);
#pragma omp barrier
		}//while over
	}//pragma over

	return counter;
}


void two_opt_random_swap(nd* min_circuit, nd cities, nd k) {
	for(nd iter=0; iter<k; iter++) {
		nd i = rand()%cities;
		nd j = rand()%cities;
		printf("%ld %ld\n",i,j);

		for(nd z=0; z<=(j-i-1)/2; z++) {
			nd temp = min_circuit[i+1+z];
			min_circuit[i+1+z] = min_circuit[j-z];
			min_circuit[j-z] = temp;
		}
	}
}
