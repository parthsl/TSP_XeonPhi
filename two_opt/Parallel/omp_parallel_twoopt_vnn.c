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
 * @Author : Parth Shah
 * TSP approximation algorithm using 2-opt method.
 * Parallelized both VNN and 2-opt technique with OpenMP.
 * 2-opt method parallelism works on chunk of some cities defined as block size.
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
	for(int i=0;i<cities-1;i++){
		int i_city = min_circuit[i];
		for(int j=i+2;j<cities;j++){
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
	omp_set_num_threads(7);

	double stime = omp_get_wtime();
	double min_tour = DBL_MAX;
	
	int* min_circuit;
	double min_tour_length=DBL_MAX;
	double tour_length;
#pragma omp parallel
{
	int id = omp_get_thread_num();
	static int block_size;
	if(id==0){
		block_size = total_cities/omp_get_num_threads();
	}
	#pragma omp barrier
	int* inner_min_circuit;
	double inner_min_tour_length = DBL_MAX;

	for(int i=block_size*id;i<((id+1)*block_size);i++){
		int* min_circuit_temp = VNN(G,total_cities,i);
		double m = find_tour_length(G,min_circuit_temp,total_cities);
		if(m < inner_min_tour_length){
			inner_min_tour_length = m;
			inner_min_circuit = min_circuit_temp;
		}
		else free(min_circuit_temp);
	}
	#pragma omp single
	{
		int rem_block_size = total_cities%omp_get_num_threads();//remaining block size
		for(int i=0;i<rem_block_size;i++){
			int* min_circuit_temp = VNN(G,total_cities,total_cities-1-i);
			double m = find_tour_length(G,min_circuit_temp,total_cities);
			if(m < inner_min_tour_length){
				inner_min_tour_length = m;
				inner_min_circuit = min_circuit_temp;
			}
			else free(min_circuit_temp);

		}
	}

	#pragma omp critical
	{
		if(min_tour_length > inner_min_tour_length){
			min_tour_length = inner_min_tour_length;
			min_circuit = inner_min_circuit;
		}
	}
}

	printf("Tour length after VNN : %lf\n",find_tour_length(G,min_circuit,total_cities));
//---------------------------------------------------------------------------------------------------
	while(1){
		tour_length = find_tour_length(G,min_circuit,total_cities);
	
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		static int block_size;
		if(id==0)
			block_size = (total_cities)/omp_get_num_threads();
		#pragma omp barrier
	
		for(int i=0;i<block_size;i++){
			if(id==omp_get_num_threads()-1){
				if(((id+1)*block_size+i)<total_cities){
					two_opt(G, min_circuit+id*block_size+i, block_size);
				}
			}
			
			else{			
				two_opt(G, min_circuit+id*block_size+i, block_size);
			}

			#pragma omp barrier //synchronize blocks movement
		}
	}
		double new_tour_length = find_tour_length(G,min_circuit,total_cities);
		if(new_tour_length == tour_length)break;
		tour_length = new_tour_length;
	}

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
