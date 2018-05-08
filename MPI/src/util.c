#include "hill_climb.h"

void donot_optimize() {
        printf("b\b");
}

nd readGraph(char* filename,nd* cities) {
        const nd BUFFER_LENGTH = 10240;
        FILE* fp = fopen(filename, "r");
        char line[BUFFER_LENGTH];
        nd graph_top;

        for(nd i=0; i<6; i++) {
                if(!fgets(line,BUFFER_LENGTH,fp))perror("Error reading file\n");
                if(strstr(line,"DIMENSION")!=NULL)
                        sscanf(line,"%*[^0123456789]%ld",cities);
        }
	
        if((*cities)!=0) {
                G=(double**)malloc(2*sizeof(double*));
		for(nd i=0;i<2;i++)
			G[i] = (double*)malloc(sizeof(double)*(*cities));
                graph_top=0;
        }
        else return EAGAIN;

        while(fgets(line,10240,fp)) {
                double x,y,t;
                if(strstr(line,"EOF")!=NULL)break;
                sscanf(line,"%lf %lf %lf",&t,&x,&y);
                G[0][graph_top] = x;
		G[1][graph_top] = y;
                graph_top++;
        }
        fclose(fp);
        return 0;
}

double find_tour_length(double** G , nd* tour, nd cities) {
        double tour_length = 0;
        for(nd i=0; i<cities-1; i++)
                tour_length += euclidean_dist(G,tour[i],tour[i+1]);

        tour_length +=euclidean_dist(G, tour[cities-1], tour[0]);

        return tour_length;
}

void print_tour(double** G, nd* min_circuit, nd total_cities) {
        for(nd i=0; i<total_cities; i++)printf("%ld ",min_circuit[i]);
        printf("\n");
        printf("Tour length after 2-opt = %lf\n",find_tour_length(G,min_circuit,total_cities));
}

void rev_arr(nd* min_circuit,nd s, nd e) {
        nd j= (e-s)/2;
#pragma omp parallel for
        for(nd i=0; i<=j; i++) {
                printf("%d\n",omp_get_num_threads());
                nd temp = min_circuit[s+i];
                min_circuit[s+i] = min_circuit[e-i];
                min_circuit[e-i] = temp;
        }
}

