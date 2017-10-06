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
 * @Author : Dhaval Shah
 * 2opt method using OpenMP
 * Each thread computes 2opt in parallel and swap egde giving best benefit
 * VNN is un-parallelized
 */
double** G;

double euclidean_dist(double x1, double y1, double x2, double y2){
        return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

typedef struct pair{
	double first,second;
}pair;



int readGraph(char* filename,int* cities){
	const int BUFFER_LENGTH = 1024;
	FILE* fp = fopen(filename, "r");
	char line[BUFFER_LENGTH];
	pair* graph;
	int graph_top;
	int graph_limit;

	for(int i=0;i<6;i++){
		fgets(line,BUFFER_LENGTH,fp);
		if(strstr(line,"DIMENSION")!=NULL)
			sscanf(line,"%*[^0123456789]%d",cities);
	}

	if(cities!=0){
		graph=(pair*)calloc(*cities,sizeof(pair));
		graph_top=0;
		graph_limit = *cities;
	}
	else return EAGAIN;

	while(fgets(line,10240,fp)){
                double x,y,t;
                if(strstr(line,"EOF")!=NULL)break;
                sscanf(line,"%lf %lf %lf",&t,&x,&y);
		graph[graph_top] = (pair){.first=x, .second=y};
		graph_top++;
        }
	fclose(fp);
	G = (double**)malloc(*cities*sizeof(double*));
	for(int i=0;i<*cities;i++)
		G[i] = (double*)calloc(*cities,sizeof(double));

	#pragma omp parallel for		
	for(int i=0;i<*cities;i++)
                for(int j=i+1;j<*cities;j++)
                        G[i][j] = G[j][i] = euclidean_dist(graph[i].first,graph[i].second, graph[j].first, graph[j].second);              
	
	free(graph);
	return 0;
}

double find_tour_length(double** G, int* tour, int cities){
	double tour_length = 0;
	for(int i=0;i<cities-1;i++)
		tour_length += G[tour[i]][tour[i+1]];

	tour_length += G[tour[cities-1]][tour[0]];

	return tour_length;
}

int* VNN(double** G, int cities, int start){
	bool* visited = (bool*)malloc(cities*sizeof(bool));
	memset(visited, 0, sizeof(bool)*cities);
	int * min_circuit = (int*) malloc(sizeof(int)*(cities+1));
	int top = 0;

	visited[start] = true;
	min_circuit[top++] = start;
	for(int i=0;i<cities-1;i++){
		double min_dist = DBL_MAX;
		int min_dist_node;
		for(int j=0;j<cities;j++){
			if(visited[j])continue;
			if(G[start][j]<min_dist){min_dist = G[start][j];min_dist_node = j;}
		}
		start = min_dist_node;
		visited[start] = true;
		min_circuit[top++] = start;
	}
	
	min_circuit[cities] = min_circuit[0];
	return min_circuit;
}

void rev_arr(int* min_circuit,int s, int e){
	int j= e;
	for(int i=s;i<j;i++){
		min_circuit[i] = min_circuit[j] + min_circuit[i] - (min_circuit[j] = min_circuit[i]);
		j--;
	}
}

void print_tour(double**G, int* min_circuit, int total_cities){
	for(int i=0;i<total_cities;i++)printf("%d ",min_circuit[i]);
	printf("\n");
	printf("Tour length after 2-opt = %lf\n",find_tour_length(G,min_circuit,total_cities));
}

void two_opt(double** G, int* min_circuit, int cities){
	double max_change = 0;
	int ic,jc;
	int counter = 0;
#pragma omp parallel 
{
	int id = omp_get_thread_num();
	int ic_local, jc_local;
	static int total_threads;
	static bool loop= true;
	static int block_size;
	#pragma omp single
	{
		total_threads = omp_get_num_threads();
		if(total_threads==0)
			block_size = 1;
		else block_size = cities/total_threads;
	}
	#pragma omp barrier
	int iend;
	if(id==total_threads-1)
		iend = cities;
	else iend = block_size*(id+1);
	#pragma omp barrier
	while(loop){	
		int i = block_size*id;
		double max_change_local = 0;
		for(;i<iend;i++){
			int i_city = min_circuit[i];
			for(int j=i+2;j<cities;j++){
				int j_city = min_circuit[j];
				int j_next_city = min_circuit[j+1];
				int i_next_city = min_circuit[i+1];
				double s_dist = G[i_city][j_city] + G[i_next_city][j_next_city];
				double f_dist = G[i_city][i_next_city] + G[j_city][j_next_city];
				if(f_dist>s_dist){
					if(f_dist-s_dist > max_change_local){
						max_change_local = f_dist-s_dist;
						ic_local = i;jc_local=j;
					}
				}
			}
		}
		#pragma omp critical
		{
			if(max_change_local > max_change){
				max_change = max_change_local;
				ic = ic_local; jc = jc_local;
			}
		}
		#pragma omp barrier
		{
			#pragma omp single
			{
				if(max_change>0){
					rev_arr(min_circuit,ic+1,jc);
					max_change = 0;
					counter++;
				}
				else{
					loop=false;
				}
			}
		}
		#pragma omp barrier
	}
}//pragma over

printf("Hill climbed = %d\n",counter);
	
}

void debug(char* filename, int cities);

int main(int argc, char** argv)
{
	int total_cities;
	srand(0);
	if(argc<2){fprintf(stderr,"Provide input file\n");return 0;}
        if(readGraph(argv[1],&total_cities)==EAGAIN)return EINVAL;

	double stime = omp_get_wtime();
	double min_tour = DBL_MAX;
	
	int* min_circuit;
	double min_tour_length=DBL_MAX;
	double tour_length;
	min_circuit = VNN(G,total_cities,0);


	printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));
//---------------------------------------------------------------------------------------------------
	tour_length = find_tour_length(G,min_circuit,total_cities);

	two_opt(G, min_circuit,total_cities);

	tour_length = find_tour_length(G,min_circuit,total_cities);

	if(tour_length<min_tour)min_tour = tour_length;

	printf("Min distance = %lf\n",min_tour);
	printf("Time taken = %lf\n",omp_get_wtime()-stime);
	free(min_circuit);
	return 0;	
}

void debug(char *filename, int cities){
	FILE * fp = fopen(filename,"r");
	int* c = (int*)malloc(cities*sizeof(int));
	for(int i=0;i<cities;i++)
	fscanf(fp,"%d",&c[i]);

	printf("fdfdf=%lf",find_tour_length(G,c,cities));
}
