#include "hill_climb.h"

nd two_opt_inline_swap(struct coords* G, nd* min_circuit, nd cities) {
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

					double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
					double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
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


nd two_opt_max_swap(struct coords* G, nd* min_circuit, nd cities) {
	double max_change = 0;
	nd ic,jc;
	nd counter = 0;
	int TILE_WIDTH = 1024;
	#pragma omp parallel
	{
		nd id = omp_get_thread_num();
		nd ic_local, jc_local;
		nd total_threads;
		static bool loop= true;

		total_threads = omp_get_num_threads();

		nd iend;
		if(id==total_threads-1)
			iend = cities;
		else iend = block_size*(id+1);
	
                nd alpha = (1-(float)id/total_threads)*(cities-2)*(cities-3);
                nd i = cities - ceil(( sqrt((alpha-6)*4 + 25)+5 )/2);
                if(id==0)i = 0;

                alpha =(1- (float)(id+1)/total_threads)*(cities-2)*(cities-3);
                iend = cities - ceil(( sqrt((alpha-6)*4 + 25) +5) /2);
                if(id+1 == total_threads)iend=cities-1;

		//To enter in while loop simultaniously
		#pragma omp barrier

		while(loop) {
	                alpha = (1-(float)id/total_threads)*(cities-2)*(cities-3);
        	        i = cities - ceil(( sqrt((alpha-6)*4 + 25)+5 )/2);
                	if(id==0)i = 0;

			double max_change_local = 0;
			for(; i<(iend-TILE_WIDTH); i+=TILE_WIDTH){
				nd j=0;
				for(j=i+2; j<(cities-2-TILE_WIDTH); j+=TILE_WIDTH){
					for(nd z=i;z<(i+TILE_WIDTH); z++){
						nd i_city = min_circuit[z];
						for(nd g=z+2;g<(j+TILE_WIDTH);g++){
							nd j_city = min_circuit[g];
		                                        nd j_next_city = min_circuit[g+1];
		                                        nd i_next_city = min_circuit[z+1];
                		                        double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
                                		        double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
		                                        if(f_dist>s_dist) {
                		                                if(f_dist-s_dist > max_change_local) {
                                		                        max_change_local = f_dist-s_dist;
                                                		        ic_local = z;
		                                                        jc_local= g;
                		                                }
                                		        }

						}
					}
				}
				for(nd z=i;z<(i+TILE_WIDTH); z++){//same i as above loops
					nd i_city = min_circuit[z];
					for(nd g=z+2;g<(cities-2);g++){//remaining j loops
						nd j_city = min_circuit[g];
                                                nd j_next_city = min_circuit[g+1];
                                                nd i_next_city = min_circuit[z+1];
                                                double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
                                                double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
                                                if(f_dist>s_dist) {
                                                        if(f_dist-s_dist > max_change_local) {
                                                                max_change_local = f_dist-s_dist;
                                                                ic_local = z;
                                                                jc_local=g;
                                                        }
                                                }
                                        }
				}
			}
			for(; i<iend;i++){//remaining i loops
				nd i_city = min_circuit[i];
				for(nd j=i+2;j<cities-2;j++){
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
