#include "hill_climb.h"

typedef long int nd;

void bruteforce(struct coords* G, nd* min_circuit, nd cities){
	nd* temp = (nd*)malloc(sizeof(nd)*cities);
	double min_dist;

	memcpy(temp,min_circuit,sizeof(nd)*cities);

	for(nd iter=0; iter<cities-1;iter++) min_dist += euclidean_dist(G[temp[iter]],G[temp[iter+1]]);

//	min_dist += euclidean_dist(G[temp[0]],G[temp[cities-1]]);//for circuit

	double dist = 0;
	
	nd fc = 1;
	for(nd i=cities;i>1;i--)fc *= i;
	nd j=1;
	nd m=0;
	
	for(nd perm=0; perm<fc;){
		nd k=0;
		while(k != fc/cities){
			while(j!= cities-1){
				dist = 0;
				for(nd iter=0; iter<cities-1;iter++) dist += euclidean_dist(G[temp[iter]],G[temp[iter+1]]);
				//dist += euclidean_dist(G[temp[cities-1]],G[temp[0]]);//fir circuit

				if(dist<min_dist){
					min_dist = dist;
					memcpy(min_circuit,temp, sizeof(nd)*cities);
				}
				
				//for(int z=0; z<cities; z++)printf("%ld ", temp[z]); //for debug
				//printf("\n");//for debug
				nd t = temp[j]; temp[j] = temp[j+1]; temp[j+1] = t;

				k++; perm++; j++;
			}
			j = 1;
		}
		m++;
		if(m==cities)break;
		k = temp[0]; temp[0] = temp[m]; temp[m] = k;
	}

	free(temp);
}


void chunk_bruteforce(struct coords* G, nd* min_circuit, nd cities, nd chunk_size){
	#pragma omp parallel
        {
                nd id = omp_get_thread_num();
                nd total_threads;
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

                nd i = block_size*id;

                for(; i<iend; i+=chunk_size){
			if(i+chunk_size < iend)
			bruteforce(G,min_circuit+i,chunk_size);
		}
	}

}
