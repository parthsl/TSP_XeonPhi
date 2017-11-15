#include "hill_climb.h"

nd two_opt_inline_swap(struct coords* G, nd* min_circuit, nd cities) {
	nd counter = 0;
	noopt* lock;
	#pragma omp parallel
	{
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
			iend = cities;
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
				double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
				double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
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
		if(id==total_threads-1)lock[id]=cities;

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

					double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
					double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
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

		#pragma omp barrier
	}//pragma over

	return counter;
}


nd two_opt_max_swap(struct coords* G, nd* min_circuit, nd cities) {
	double max_change = 0;
	nd ic,jc;
	nd counter = 0;
	#pragma omp parallel
	{
		nd id = omp_get_thread_num();
		nd ic_local, jc_local;
		nd total_threads;
		static bool loop= true;
		nd block_size;

		total_threads = omp_get_num_threads();
		if(total_threads==0)
			block_size = 1;
		else block_size = cities/total_threads;

		nd iend;
		if(id==total_threads-1)
			iend = cities;
		else iend = block_size*(id+1);

		#pragma omp barrier

		while(loop) {
			nd i = block_size*id;
			double max_change_local = 0;
			for(; i<iend; i++) {
				nd i_city = min_circuit[i];
				for(nd j=i+2; j<cities-1; j++) {
					nd j_city = min_circuit[j];
					nd j_next_city = min_circuit[j+1];
					nd i_next_city = min_circuit[i+1];
					double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
					double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
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

			static nd j;
			#pragma omp single
			j = (jc-ic-1)/2;

			#pragma omp barrier
			if(max_change>0) {
				#pragma omp for
				for(i=0; i<=j; i++) {
					nd temp = min_circuit[ic+1+i];
					min_circuit[ic+1+i] = min_circuit[jc-i];
					min_circuit[jc-i] = temp;
				}
			}
			else loop=false;

			#pragma omp barrier
			max_change = 0;
			#pragma omp single
			counter++;
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
