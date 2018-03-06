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


nd two_opt_max_swap(struct coords* G, nd* min_circuit, const nd cities) {
	double max_change = 0;
	nd ic,jc;
	nd counter = 0;
	bool loop= true;
	double precal_distance[cities];

	#pragma omp parallel
	{
		nd id = omp_get_thread_num();
		nd ic_local=0, jc_local=0;
		nd total_threads;
		nd i,iend;
		nd alpha;
		const long long int VS = 64;
		int VVS = VS;
		double avx_ed[VS];
		double avx_pre[VS];
#if defined(__ibmxl__) || defined(__powerpc__)
		vector double x = {0,0},y={0,0};
#elif defined(__INTEL_COMPILER)
		double x[2]={0,0}, y[2]={0,0};
#endif

		total_threads = omp_get_num_threads();
		alpha =(1- (float)(id+1)/total_threads)*(cities-2)*(cities-3);
		iend = cities - ceil(( sqrt((alpha-6)*4 + 25) +5) /2);
		if(id+1 == total_threads)iend=cities-1;

		#pragma omp for schedule(static)
		for(i=0; i<cities; i++) precal_distance[i] = euclidean_dist(G[min_circuit[i]],G[min_circuit[i+1]]);

		//To enter in while loop simultaniously
		#pragma omp barrier

		alpha = (1-(float)id/total_threads)*(cities-2)*(cities-3);
		alpha = cities - ceil(( sqrt((alpha-6)*4 + 25)+5 )/2);
		if(id==0)alpha=0;
		while(loop) {
			i = alpha;
			double max_change_local = 0;

			for(; i<iend; i++) {
				nd i_city = min_circuit[i];
				nd i_next_city = min_circuit[i+1];
				nd j = i+2;

				for(VVS=VS; VVS>=8; VVS=VVS/2) {
					for(; j<(cities-1-VVS); j+=VVS) {
#if defined(__ibmxl__) || defined(__powerpc__)
						for(nd jj=0; jj<VVS; jj++) {
							avx_ed[jj] = squared_dist(G[i_city],G[min_circuit[j+jj]]);
							avx_pre[jj] = squared_dist(G[i_next_city],G[min_circuit[j+jj+1]]);
						}
						vsqrt(avx_ed,avx_ed,&VVS);
						vsqrt(avx_pre,avx_pre,&VVS);
						for(nd jj=0; jj<VVS; jj++) avx_ed[jj] = avx_ed[jj] + avx_pre[jj];
#elif defined(__INTEL_COMPILER)
						#pragma omp simd
						for(nd jj=0; jj<VVS; jj++) {
							avx_ed[jj] = squared_dist(G[i_city],G[min_circuit[j+jj]]);
							avx_pre[jj] = squared_dist(G[i_next_city],G[min_circuit[j+jj+1]]);
						}
						vdSqrt(VVS,avx_ed,avx_ed);
						vdSqrt(VVS,avx_pre,avx_pre);
						#pragma omp simd
						for(nd jj=0; jj<VVS; jj++) avx_ed[jj] = avx_ed[jj] + avx_pre[jj];
#else
						for(nd jj=0; jj<VVS; jj++) {
							avx_ed[jj] = euclidean_dist(G[i_city],G[min_circuit[j+jj]]) + euclidean_dist(G[i_next_city],G[min_circuit[j+jj+1]]);
						}
#endif

						#pragma omp simd
						for(nd jj=0; jj<VVS; jj++) {
							avx_pre[jj] = precal_distance[i] + precal_distance[j+jj];
						}
						#pragma omp simd
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

				//Peeled loop.
				for(; j<cities-1; j++) {
					nd j_city = min_circuit[j];
					nd j_next_city = min_circuit[j+1];
#if defined(__ibmxl__) || defined(__powerpc__)
					x[0] = squared_dist(G[i_city],G[j_city]);
					x[1] = squared_dist(G[i_next_city],G[j_next_city]);
					y = sqrtd2(x);
					double s_dist = y[0] + y[1];
#elif defined(__INTEL_COMPILER)
					x[0] = squared_dist(G[i_city],G[j_city]);
					x[1] = squared_dist(G[i_next_city],G[j_next_city]);
					vdSqrt(2,x,y);
					double s_dist = y[0] + y[1];
#else
					double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
#endif
					double f_dist = precal_distance[i]+precal_distance[j];//euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
					if(f_dist>s_dist) {
						if(f_dist-s_dist > max_change_local) {
							max_change_local = f_dist-s_dist;
							ic_local = i;
							jc_local = j;
						}
					}
				}
			}

			static nd j;
			#pragma omp critical
			{
				if(max_change_local > max_change) {
					max_change = max_change_local;
					ic = ic_local;
					jc = jc_local;
					j = (jc-ic-1)/2;
				}
			}

			#pragma omp barrier
			if(max_change>0) {
				#pragma omp for schedule(static)
				for(i=0; i<=j; i++) {
					nd temp = min_circuit[ic+1+i];
					min_circuit[ic+1+i] = min_circuit[jc-i];
					min_circuit[jc-i] = temp;
				}
				#pragma omp for schedule(static)
				for(i=0; i<=j; i++) {
					nd temp = precal_distance[ic+i];
					precal_distance[ic+i] = precal_distance[jc-i];
					precal_distance[jc-i] = temp;
				}
#if defined(__ibmxl__) || defined(__powerpc__)
				x[0] = squared_dist(G[min_circuit[ic]],G[min_circuit[ic+1]]);
				x[1] = squared_dist(G[min_circuit[jc]], G[min_circuit[jc+1]]);
				y = sqrtd2(x);
				precal_distance[ic] = y[0];
				precal_distance[jc] = y[1];
#elif defined(__INTEL_COMPILER)
				x[0] = squared_dist(G[min_circuit[ic]],G[min_circuit[ic+1]]);
				x[1] = squared_dist(G[min_circuit[jc]], G[min_circuit[jc+1]]);
				vdSqrt(2,x,y);
				precal_distance[ic] = y[0];
				precal_distance[jc] = y[1];

#else
				precal_distance[ic] = euclidean_dist(G[min_circuit[ic]],G[min_circuit[ic+1]]);
				precal_distance[jc] = euclidean_dist(G[min_circuit[jc]],G[min_circuit[jc+1]]);
#endif
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
		printf("%ld%ld\n",i,j);

		for(nd z=0; z<=(j-i-1)/2; z++) {
			nd temp = min_circuit[i+1+z];
			min_circuit[i+1+z] = min_circuit[j-z];
			min_circuit[j-z] = temp;
		}
	}
}


/*
 * Single Thread Execution of  two_opt_max_swap() function
 * Optimised for XLC compiler using MASS library for sqrt vector instructions
 */
nd two_opt_max_swap_single(struct coords* G, nd* min_circuit, nd cities) {
	double max_change = 0;
	nd counter = 0;
	bool loop = true;
	nd ic = 0,jc = 0;
	const long long int VS = 64;//Vector Size: Length of vectorised operation queue.
	int VVS = (int)cities;//Variable Vector Size: Length of vectorising operation queue varies with each iteration.
	double avx_ed[VS];
	double avx_pre[VS];
#if defined(__ibmxl__) || defined(__powerpc__)
	vector double x={0,0},y={0,0};
#elif defined(__INTEL_COMPILER)
	double x[2]={0,0}, y[2]={0,0};
	printf("Detected ICC compiled binary\n");
#endif

	double precal_distance[cities];

#if defined(__ibmxl__) || defined(__powerpc__)
	for(nd i=0; i<cities-1; i++) {
		precal_distance[i] = squared_dist(G[min_circuit[i]],G[min_circuit[i+1]]);
	}
	vsqrt(precal_distance, precal_distance, &VVS);
#else
	for(nd i=0; i<cities-1; i++) {
		precal_distance[i] = euclidean_dist(G[min_circuit[i]],G[min_circuit[i+1]]);
	}
#endif

	while(loop) {

		nd i=0,j = 0;
		for(; i<cities-2; i++) {
			VVS = VS;

			nd i_city = min_circuit[i];
			nd i_next_city = min_circuit[i+1];
			j = i+2;

			/*Vectorising operation to find maximum swap
			 * Reducing Vector Size(VVS) to calculate for peeled loop in each iteration
			 * VVS reduced till 2 due to 128bit vector length in POWER
			 * Change lower limit to 4/8 for intel machines.
			 */
			for(VVS=VS; VVS>=2; VVS/=2) {

				for(; j<(cities-1-VVS); j+=VVS) {
#if defined(__ibmxl__) || defined(__powerpc__)
					for(nd jj=0; jj<VVS; jj++) {
						avx_ed[jj] = squared_dist(G[i_city],G[min_circuit[j+jj]]);
						avx_pre[jj] = squared_dist(G[i_next_city],G[min_circuit[j+jj+1]]);
					}
					vsqrt(avx_ed,avx_ed,&VVS);
					vsqrt(avx_pre,avx_pre,&VVS);
					for(nd jj=0; jj<VVS; jj++) avx_ed[jj] = avx_ed[jj] + avx_pre[jj];
#elif defined(__INTEL_COMPILER)
					for(nd jj=0; jj<VVS; jj++) {
						avx_ed[jj] = squared_dist(G[i_city],G[min_circuit[j+jj]]);
						avx_pre[jj] = squared_dist(G[i_next_city],G[min_circuit[j+jj+1]]);
					}
					vdSqrt(VVS,avx_ed,avx_ed);
					vdSqrt(VVS,avx_pre,avx_pre);
					for(nd jj=0; jj<VVS; jj++) avx_ed[jj] = avx_ed[jj] + avx_pre[jj];
#else
					for(nd jj=0; jj<VVS; jj++) {
						avx_ed[jj] = euclidean_dist(G[i_city],G[min_circuit[j+jj]]) + euclidean_dist(G[i_next_city],G[min_circuit[j+jj+1]]);
					}
#endif

					for(nd jj=0; jj<VVS; jj++) {
						avx_pre[jj] = precal_distance[i] + precal_distance[j+jj];
					}
					for(nd jj=0; jj<VVS; jj++) {
						avx_pre[jj] = avx_pre[jj] - avx_ed[jj];
					}
					for(nd jj=0; jj<VVS; jj++) {
						if(avx_pre[jj] > max_change) {
							max_change = avx_pre[jj];
							ic = i;
							jc = j+jj;
						}
					}
				}
			}

			//Peeled loop.
			for(; j<cities-1; j++) {
				nd j_city = min_circuit[j];
				nd j_next_city = min_circuit[j+1];
#if defined(__ibmxl__) || defined(__powerpc__)
				x[0] = squared_dist(G[i_city],G[j_city]);
				x[1] = squared_dist(G[i_next_city],G[j_next_city]);
				y = sqrtd2(x);
				double s_dist = y[0] + y[1];
#elif defined(__INTEL_COMPILER)
				x[0] = squared_dist(G[i_city],G[j_city]);
				x[1] = squared_dist(G[i_next_city],G[j_next_city]);
				vdSqrt(2,x,y);
				double s_dist = y[0] + y[1];
#else
				double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
#endif
				double f_dist = precal_distance[i]+precal_distance[j];//euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
				if(f_dist>s_dist) {
					if(f_dist-s_dist > max_change) {
						max_change = f_dist-s_dist;
						ic = i;
						jc = j;
					}
				}
			}
		}

		j = (jc-ic-1)/2;


		if(max_change>0) {
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
			x[0] = squared_dist(G[min_circuit[ic]],G[min_circuit[ic+1]]);
			x[1] =  squared_dist(G[min_circuit[jc]], G[min_circuit[jc+1]]);
			y = sqrtd2(x);
			precal_distance[ic] = y[0];
			precal_distance[jc] = y[1];
#elif defined(__INTEL_COMPILER)
			x[0] = squared_dist(G[min_circuit[ic]],G[min_circuit[ic+1]]);
			x[1] = squared_dist(G[min_circuit[jc]], G[min_circuit[jc+1]]);
			vdSqrt(2,x,y);
			precal_distance[ic] = y[0];
			precal_distance[jc] = y[1];
#else
			precal_distance[ic] = euclidean_dist(G[min_circuit[ic]],G[min_circuit[ic+1]]);
			precal_distance[jc] = euclidean_dist(G[min_circuit[jc]],G[min_circuit[jc+1]]);
#endif
		}
		else loop=false;

		max_change = 0;
		counter++;
	}//while over

	return counter;
}

