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
	nd sub_block_size = cities/(nd)num_procs;
	nd sub_block_start = num_local*sub_block_size;
	nd sub_block_end = sub_block_size*(num_local+1);
	if(num_local == num_procs-1){
		sub_block_end = cities;
		sub_block_size = sub_block_end - sub_block_start;
	}

#ifdef DEBUG
	printf("%d: block_size=%ld\n",num_local,sub_block_size);
#endif

	//for Openmp
	double max_change = 0;
	nd ic=0,jc=0;
	nd counter = 0;

#pragma omp parallel
	{
		nd id = omp_get_thread_num();
		nd total_threads = omp_get_num_threads();
		nd iend;
		nd tile_size;

		enum loopFlags{ Proceed, Abort};
		static int loop_flag= Proceed;

		if(total_threads==0)
			tile_size = 1;
		else tile_size = sub_block_size/total_threads;

		if(id==total_threads-1)
			iend = sub_block_end;

		else iend = sub_block_start + tile_size*(id+1);

#ifdef DEBUG
		printf("%d:%ld block_start=%ld end=%ld\n",num_local,id,sub_block_start+tile_size*id, iend);
#endif

#pragma omp barrier

		while(loop_flag != Abort) {
			nd ic_local=0, jc_local=0;
			nd i = sub_block_start + tile_size*id;
			double max_change_local = 0;
			for(; i<iend; i++) { //goto iend as OutOfRange Checking is done by inner j loop
				nd i_city = min_circuit[i];
				for(nd j=i+2; j<cities-1; j++) {
					nd j_city = min_circuit[j];
					nd j_next_city = min_circuit[j+1];
					nd i_next_city = min_circuit[i+1];
					if(i>=cities || j>=cities)fprintf(stderr,"Invalid i nad j values\n");
					double s_dist = euclidean_dist(G, i_city, j_city) + euclidean_dist(G,i_next_city, j_next_city);
					double f_dist = euclidean_dist(G, i_city, i_next_city) + euclidean_dist(G, j_city, j_next_city);
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

			/* Print error on invalid city coordinates obtained*/
			if(ic >= cities || jc >= cities)
				fprintf(stderr, "Invalid values%ld %ld %lf\n",ic,jc,max_change);

#pragma omp barrier
			if(id==0)
				MPI_Barrier(MPI_COMM_WORLD);
#pragma omp barrier

			/*----------------Find best swap among all nodes--------*/
			if(id==0){
				double* max_change_procs;
				nd selected_coords[2];
				nd selected_nodes_swap = 0;
				if(num_local == master) //master will save all distance computed by other nodes
					max_change_procs = (double*)malloc(num_procs*sizeof(double));

				MPI_Gather(&max_change, 1, MPI_DOUBLE, max_change_procs, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);//to get all Max Swaps


				if(num_local == master){ //find node with maxswap
					for(nd i = 0; i<num_procs-1; i++){
						if(max_change < max_change_procs[i]){
							selected_nodes_swap = i;
							max_change = max_change_procs[i];
						}
					}

					#ifdef DEBUG
					printf("Max Change selected node = %ld with distance=%lf\n",selected_nodes_swap,max_change);
					#endif
					free(max_change_procs);
				}

				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Bcast(&max_change, 1, MPI_DOUBLE, master, MPI_COMM_WORLD); //Send Max Change Distance
				MPI_Bcast(&selected_nodes_swap, 1, MPI_LONG_LONG, master, MPI_COMM_WORLD); //Send Selected node in MPI Cluster giving best swap

				selected_coords[0] = ic;
				selected_coords[1] = jc;
				if(num_local == selected_nodes_swap && num_local != master ){ // Send ic and jc from Selected node to MASTER node
#ifdef DEBUG
					printf("%d sends ic=%ld jd=%ld to master\n",num_local,selected_coords[0],selected_coords[1]);
#endif
					MPI_Send(selected_coords,2,MPI_LONG_LONG, master, 0, MPI_COMM_WORLD);
				}
				else if(num_local != selected_nodes_swap && num_local==master){ //Master receives the selected ic and jc cities
					MPI_Recv(selected_coords, 2, MPI_LONG_LONG, selected_nodes_swap,0,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef DEBUG
					printf("Master recieves ic=%ld jc=%ld\n",selected_coords[0],selected_coords[1]);
#endif
				}

#ifdef DEBUG
				if(num_local==master)
					printf("Broadcasting selected ic and jc to all other nodes\n");
#endif

				MPI_Bcast(selected_coords, 2, MPI_LONG_LONG, master, MPI_COMM_WORLD); //Broadcast Selected ic and jc to each nodes
				ic = selected_coords[0];
				jc = selected_coords[1];
				if(num_local == master && max_change > 0)
					loop_flag = Proceed;
				else loop_flag = Abort;

#ifdef DEBUG
				if(num_local==master)
					printf("Broadcasting loop flag=%d\n",loop_flag);
#endif

				MPI_Bcast(&loop_flag, 1, MPI_INT, master, MPI_COMM_WORLD); // Send whether to proceed or stop

				/*------------------swap array---------------------------*/

				if(loop_flag == Proceed){
					static nd j;
					j = (jc-ic-1)/2;
					for(i=0; i<=j; i++) {
						nd temp = min_circuit[ic+1+i];
						min_circuit[ic+1+i] = min_circuit[jc-i];
						min_circuit[jc-i] = temp;
					}
					max_change = 0;
					ic=0;jc=0;
					counter++;
				}
				else{
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

/*-------------MPICH------------------version for Max_swap------------------------*/
nd two_opt_tiled_max_swap(double** G, nd* min_circuit, nd cities, int num_procs, int num_local) {
        //for Mpich
        nd sub_block_size = cities/(nd)num_procs;
        nd sub_block_start = num_local*sub_block_size;
        nd sub_block_end = sub_block_size*(num_local+1);
        if(num_local == num_procs-1){
                sub_block_end = cities;
                sub_block_size = sub_block_end - sub_block_start;
        }

#ifndef DEBUG
        printf("%d: block_size=%ld\n",num_local,sub_block_size);
#endif

        //for Openmp
        double max_change = 0;
        nd ic=0,jc=0;
        nd counter = 0;

#pragma omp parallel
        {
                nd id = omp_get_thread_num();
                nd total_threads = omp_get_num_threads();
                nd iend;
                nd tile_size;

                enum loopFlags{ Proceed, Abort};
                static int loop_flag= Proceed;

                if(total_threads==0)
                        tile_size = 1;
                else tile_size = sub_block_size/total_threads;

                if(id==total_threads-1)
                        iend = sub_block_end;

                else iend = sub_block_start + tile_size*(id+1);

#ifdef DEBUG
                printf("%d:%ld block_start=%ld end=%ld\n",num_local,id,sub_block_start+tile_size*id, iend);
#endif

#pragma omp barrier
                while(loop_flag != Abort){
                        nd ic_local=0, jc_local=0;
                        nd i = sub_block_start + tile_size*id;
                        double max_change_local = 0;
                        for(; i<iend; i++) { //goto iend as OutOfRange Checking is done by inner j loop
                                nd i_city = min_circuit[i];
                                for(nd j=i+2; j<sub_block_end-1; j++) {
                                        nd j_city = min_circuit[j];
                                        nd j_next_city = min_circuit[j+1];
                                        nd i_next_city = min_circuit[i+1];
                                        if(i>=cities || j>=cities)fprintf(stderr,"Invalid i nad j values\n");
                                        double s_dist = euclidean_dist(G, i_city, j_city) + euclidean_dist(G,i_next_city, j_next_city);
                                        double f_dist = euclidean_dist(G, i_city, i_next_city) + euclidean_dist(G, j_city, j_next_city);
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

                        /* Print error on invalid city coordinates obtained*/
                        if(ic >= cities || jc >= cities)
                                fprintf(stderr, "Invalid values%ld %ld %lf\n",ic,jc,max_change);

#pragma omp barrier

                        /*-----------------------------------------------*/
#pragma omp single
                        if(max_change>0){
                                static nd j;
                                j = (jc-ic-1)/2;
                                for(i=0; i<=j; i++) {
                                        nd temp = min_circuit[ic+1+i];
                                        min_circuit[ic+1+i] = min_circuit[jc-i];
                                        min_circuit[jc-i] = temp;
                                }
                                max_change = 0;
                                ic=0;jc=0;
                                counter++;
                        }else{
                                loop_flag = Abort;
                        }

                }

                /*Now send calculated array to master node*/

#pragma omp barrier
                if(id==0)
                        MPI_Barrier(MPI_COMM_WORLD);
#pragma omp barrier
                /*MPI_THREAD_FUNNELED*/
                if(id==0){
                        #ifdef DEBUG
                        printf("%d: block start=%ld end=%ld\n",num_local,sub_block_start,sub_block_end);
                        #endif

                        sub_block_size = cities/num_procs;
                        nd* temp_circuit = (nd*)malloc(sub_block_size*sizeof(nd));
                        nd* recv_temp_circuit;

                        memcpy(temp_circuit, min_circuit+sub_block_start, sub_block_size*sizeof(nd));
                        if(num_local == master){
                                recv_temp_circuit = min_circuit;//(nd*)malloc(cities*sizeof(nd));
                        }

                        /*------------All gather to get Optimum Sub-set of Problem--------------------*/
                        MPI_Gather(&temp_circuit[0], sub_block_size, MPI_LONG_LONG, recv_temp_circuit, sub_block_size, MPI_LONG_LONG, master, MPI_COMM_WORLD);


                        if(num_local == num_procs-1 && cities%num_procs != 0){
                                //printf("%d: Sending remaining coords with nelements=%ld\n",num_local, cities-sub_block_size*num_procs);
                                MPI_Send(min_circuit+sub_block_start+sub_block_size, cities-sub_block_size*num_procs, MPI_LONG_LONG, master, 0, MPI_COMM_WORLD);
                        }
                        else if(num_local == master && cities%num_procs != 0){
                                //printf("%d: Receiving remaining coords%ld\n",num_local, cities-sub_block_size*num_procs);
                                MPI_Recv(recv_temp_circuit+sub_block_size*(num_procs), cities-(sub_block_size*num_procs), MPI_LONG_LONG, num_procs-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                }
        #pragma omp barrier
        }//pragma over

        return counter;
}
