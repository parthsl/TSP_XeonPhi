#include "hill_climb.h"

nd* VNNp(double** G, nd cities, nd start) {
	bool* visited = (bool*)calloc(cities,sizeof(bool));
	nd *min_circuit = (nd*)malloc(sizeof(nd)*(cities+1));
	min_circuit[0] = start;
	min_circuit[cities] = start;
	visited[start] = true;

	#pragma omp parallel
	{
		static nd i=1;
		static double max_distance = DBL_MAX;
		static nd node_selected;
		#pragma omp barrier
		while(i<cities) {
			nd id = omp_get_thread_num();
			nd threads = omp_get_num_threads();
			nd bs = cities/threads;
			nd jend = (id+1)*bs;
			if(id==threads-1)jend=cities;
			double max_dist_local = DBL_MAX;
			nd selected_node=0;

			for(nd j=id*bs; j<jend; j++) {
				if(visited[j])continue;
				double dist = euclidean_dist(G, start, j);
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

//Simple VNN un-parallelized
nd* VNN(double** G, nd cities, nd start) {
	bool* visited = (bool*)malloc(cities*sizeof(bool));
	memset(visited, 0, sizeof(bool)*cities);
	nd * min_circuit = (nd*) malloc(sizeof(nd)*(cities+1));
	nd top = 0;

	visited[start] = true;
	min_circuit[top++] = start;
	for(int i=0; i<cities-1; i++) {
		double min_dist = DBL_MAX;
		int min_dist_node = 0;
		for(int j=0; j<cities; j++) {
			if(visited[j])continue;
			if(euclidean_dist(G, start, j)<min_dist) {
				min_dist = euclidean_dist(G, start, j);
				min_dist_node = j;
			}
		}
		start = min_dist_node;
		visited[start] = true;
		min_circuit[top++] = start;
	}

	min_circuit[cities] = min_circuit[0];
	return min_circuit;
}

