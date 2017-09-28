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
 * Share distance graph among all threads.
 * Each threads works independently with using Variable Nearest Neighbour method to generate initial route.
 * Each thread will then use initial route in 2-opt method to find optimal path heuristically.
 */ 

double** G;

double euclidean_dist(double x1, double y1, double x2, double y2){
        return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

typedef struct pair{
	double first,second;
}pair;


/*
 * Reads graph from a file having coordinates of the cities.
 * Stores number of cities in the second arguments.
 */
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
		
	for(int i=0;i<*cities;i++)
                for(int j=i+1;j<*cities;j++)
                        G[i][j] = G[j][i] = euclidean_dist(graph[i].first,graph[i].second, graph[j].first, graph[j].second);              
	
	free(graph);
	return 0;
}

//Given a tour and distance matrix will give tour length.
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
	for(int i=0;i<cities-1;i++){
		int i_city = min_circuit[i];
		int z = 1;
		if(i==0)z=0;
		for(int j=i+2;j<cities-1+z;j++){
			int j_city = min_circuit[j];
			int j_next_city = min_circuit[j+1];
			int i_next_city = min_circuit[i+1];
			double s_dist = G[i_city][j_city] + G[i_next_city][j_next_city];
			double f_dist = G[i_city][i_next_city] + G[j_city][j_next_city];
			if(f_dist>s_dist){
				rev_arr(min_circuit,i+1,j);
			}
		}
	}
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
#pragma omp parallel
{
	int* min_circuit;
	double min_tour_length=DBL_MAX;
	double tour_length;
	int id = omp_get_thread_num();
	if(id<total_cities){
		min_circuit = VNN(G,total_cities,id);
		printf("Tour length after VNN(id=%d) : %lf\n",id,find_tour_length(G,min_circuit,total_cities));
		while(1){
			tour_length = find_tour_length(G,min_circuit,total_cities);
			two_opt(G,min_circuit,total_cities);
			double new_tour_length = find_tour_length(G,min_circuit,total_cities);
			if(new_tour_length == tour_length)break;
			tour_length = new_tour_length;
		}
		#pragma omp critical
		if(tour_length<min_tour)min_tour = tour_length;
	}
}
	printf("Min distance = %lf\n",min_tour);
	printf("Time taken = %lf\n",omp_get_wtime()-stime);
	return 0;	
}
