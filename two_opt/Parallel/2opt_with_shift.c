#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits.h>
#include<stdbool.h>
#include<string.h>
#include<errno.h>
#include<float.h>
#include<omp.h>

/*
 * @Author : Parth Shah<parths1229@gmail.com>A
 *
 * Parallel VNN + 2opt(w/ Tiling and inline swapping) method combined method with OpenMP giving parallel execution of both in optimal way with co-ordinates stored for each city.
 * 2-opt works in tiles of data. Each thread has their own working set and tries to swap in that workspace without waiting for other threads.
 * 2-opt uses inline swapping technique : Swap all the edges giving the benefit.
 * Results : Speedup in VNN calculation. 2opt is more faster than w/o Tiling strategy + maxSwap technique but has increase in error rate. Increase in threads can give more faster result but increase in error from optimal distance too. Able to load and process 85k cities in few seconds.
 */

struct coords {
    double x,y;
};

struct coords* G;

typedef volatile long int noopt;

#define euclidean_dist(i, j) \
	sqrt( (i.x-j.x)*(i.x-j.x) + (i.y-j.y)*(i.y-j.y))

void donot_optimize() {
    printf("b\b");
}


long int readGraph(char* filename,long int* cities) {
    const long int BUFFER_LENGTH = 10240;
    FILE* fp = fopen(filename, "r");
    char line[BUFFER_LENGTH];
    long int graph_top;

    for(long int i=0; i<6; i++) {
        if(!fgets(line,BUFFER_LENGTH,fp))perror("Error reading file\n");
        if(strstr(line,"DIMENSION")!=NULL)
            sscanf(line,"%*[^0123456789]%ld",cities);
    }

    if(cities!=0) {
        G=(struct coords*)calloc(*cities,sizeof(struct coords));
        graph_top=0;
    }
    else return EAGAIN;

    while(fgets(line,10240,fp)) {
        double x,y,t;
        if(strstr(line,"EOF")!=NULL)break;
        sscanf(line,"%lf %lf %lf",&t,&x,&y);
        G[graph_top] = (struct coords) {
            .x=x, .y=y
        };
        graph_top++;
    }
    fclose(fp);
    return 0;
}

double find_tour_length(struct coords* G , long int* tour, long int cities) {
    double tour_length = 0;
    for(long int i=0; i<cities-1; i++)
        tour_length += euclidean_dist(G[tour[i]],G[tour[i+1]]);

    tour_length +=euclidean_dist(G[tour[cities-1]],G[tour[0]]);

    return tour_length;
}


long int* VNNp(struct coords* G, long int cities, long int start) {
    bool* visited = (bool*)calloc(cities,sizeof(bool));
    long int *min_circuit = (long int*)malloc(sizeof(long int)*cities);
    min_circuit[0] = start;
    visited[start] = true;

    #pragma omp parallel
    {
        static long int i=1;
        static double max_distance = DBL_MAX;
        static long int node_selected;

        #pragma omp barrier
        while(i<cities) {
            long int id = omp_get_thread_num();
            long int threads = omp_get_num_threads();
            long int bs = cities/threads;
            long int jend = (id+1)*bs;
            if(id==threads-1)jend=cities;
            double max_dist_local = DBL_MAX;
            long int selected_node=0;
            for(long int j=id*bs; j<jend; j++) {
                if(visited[j])continue;
                double dist = euclidean_dist(G[start],G[j]);
                if(dist<max_dist_local) {
                    max_dist_local = dist;
                    selected_node = j;
                }
            }
            #pragma omp critical
            {
                if(max_dist_local < max_distance) {
                    max_distance = max_dist_local;
                    node_selected = selected_node;
                }
            }

            #pragma omp barrier
            #pragma omp single
            {
                start = node_selected;
                visited[start] = true;
                min_circuit[i] = start;
                max_distance = DBL_MAX;
                i++;
            }
            #pragma omp barrier
        }
    }
    return min_circuit;
}

long int* VNN(struct coords* G, long int cities, long int start) {
    bool* visited = (bool*)malloc(cities*sizeof(bool));
    memset(visited, 0, sizeof(bool)*cities);
    long int * min_circuit = (long int*) malloc(sizeof(long int)*cities);
    long int top = 0;

    visited[start] = true;
    min_circuit[top++] = start;
    for(long int i=0; i<cities-1; i++) {
        double min_dist = DBL_MAX;
        long int min_dist_node;
        #pragma omp parallel
        {
            long int id = omp_get_thread_num();
            long int threads = omp_get_num_threads();
            long int bs = cities/threads;
            long int jend = (id+1)*bs;
            if(id==threads-1)jend = cities;
            double min_dist_local = DBL_MAX;
            long int min_node_local=0;

            for(long int j=id*bs; j<jend; j++) {
                if(visited[j])continue;
                if(euclidean_dist(G[start],G[j])<min_dist_local) {
                    min_dist_local = euclidean_dist(G[start],G[j]);
                    min_node_local = j;
                }
            }

            #pragma omp critical
            {
                if(min_dist_local < min_dist) {
                    min_dist = min_dist_local;
                    min_dist_node = min_node_local;
                }
            }
        }
        start = min_dist_node;
        visited[start] = true;
        min_circuit[top++] = start;

    }

    return min_circuit;
}

inline void rev_arr(long int* min_circuit,long int s, long int e) {
    long int j= (e-s)/2;
    #pragma omp parallel for
    for(long int i=0; i<=j; i++) {
        printf("%d\n",omp_get_num_threads());
        long int temp = min_circuit[s+i];
        min_circuit[s+i] = min_circuit[e-i];
        min_circuit[e-i] = temp;
    }
}

void print_tour(struct coords* G, long int* min_circuit, long int total_cities) {
    for(long int i=0; i<total_cities; i++)printf("%ld ",min_circuit[i]);
    printf("\n");
    printf("Tour length after 2-opt = %lf\n",find_tour_length(G,min_circuit,total_cities));
}

void two_opt(struct coords* G, long int* min_circuit, long int cities) {
    long int counter = 0;
    noopt* lock;
    #pragma omp parallel
    {
        long int id = omp_get_thread_num();
        long int total_threads;
        long int block_size;

        total_threads = omp_get_num_threads();
        if(total_threads==0)
            block_size = 1;
        else block_size = cities/total_threads;

        #pragma omp single
        printf("Thread spawn : %ld\n",total_threads);

        #pragma omp single
        lock = (noopt*)calloc(total_threads,sizeof(noopt));

        long int iend;
        if(id==total_threads-1)
            iend = cities;
        else iend = block_size*(id+1);
        #pragma omp barrier

        long int i = block_size*id;
        double max_change_local = 0;
        for(; i<iend; i++) {
            long int i_city = min_circuit[i];
            for(long int j=i+2; j<iend; j++) {
                long int j_city = min_circuit[j];
                long int j_next_city = min_circuit[j+1];
                long int i_next_city = min_circuit[i+1];
                double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
                double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
                if(f_dist>s_dist) {
                    if(f_dist-s_dist > max_change_local) {
                        for(long int z=0; z<=(j-i-1)/2; z++) {
                            long int temp = min_circuit[i+1+z];
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

            long int j = block_size*(id+1);

            for(long int g=block_size*id+1; g<iend; g++) { // Runs through whole block

                while(lock[id+1]<j);
                long int j_next_city = min_circuit[j+1];
                for(i=g; i<iend; i++) // runs i in whole block but taking j+1 as const
                {
                    long int i_city = min_circuit[i];
                    long int i_next_city = min_circuit[i+1];
                    long int j_city = min_circuit[j];

                    double s_dist = euclidean_dist(G[i_city],G[j_city]) + euclidean_dist(G[i_next_city],G[j_next_city]);
                    double f_dist = euclidean_dist(G[i_city],G[i_next_city]) + euclidean_dist(G[j_city],G[j_next_city]);
                    if(f_dist>s_dist) {
                        if(f_dist-s_dist > max_change_local) {
                            for(long int z=0; z<=(j-i-1)/2; z++) {
                                long int temp = min_circuit[i+1+z];
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

    printf("Hill climbed = %ld\n",counter);

}


int main(int argc, char** argv) // argv1 = filename argv2 = threads argv3=debug
{
    long int total_cities;
    srand(0);
    if(argc<2) {
        fprintf(stderr,"Provide input file\n");
        return 0;
    }
    if(readGraph(argv[1],&total_cities)==EAGAIN)return EINVAL;
    if(argc>2)omp_set_num_threads(atoi(argv[2]));

    double stime = omp_get_wtime();
    double min_tour = DBL_MAX;

    long int* min_circuit;
    double tour_length;
    min_circuit = VNNp(G,total_cities,0);
    printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));
    //---------------------------------------------------------------------------------------------------
    tour_length = find_tour_length(G,min_circuit,total_cities);

    two_opt(G, min_circuit,total_cities);

    tour_length = find_tour_length(G,min_circuit,total_cities);

    if(tour_length<min_tour)min_tour = tour_length;

    printf("Min distance = %lf\n",min_tour);
    printf("Time taken = %lf\n",omp_get_wtime()-stime);
#ifdef DEBUG
    print_tour(G,min_circuit,total_cities);
#endif
    free(min_circuit);
    return 0;
}
