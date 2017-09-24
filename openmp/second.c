#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define numthreads 4


const double strips = 10000000.0;

int main(int argc, char* argv[]){
	double st = omp_get_wtime();
	double sum=0, stripwidth=1/strips;
	int threads;

	if(argc<2)
	omp_set_num_threads(numthreads);
	else 
	omp_set_num_threads(atoi(argv[1]));
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		if(id==0)threads = omp_get_num_threads();
		double summ=0;
		#pragma omp barrier
		for(int i=id;i<strips;i+=threads){
			double x = i*stripwidth;
			summ = summ + (4/(1+x*x));
		}

		#pragma omp atomic
		sum = sum+summ;
	}

	printf("Pi = %lf\n",sum/strips);

	printf("Time Taken = %g\n",omp_get_wtime()-st);
	return 0;
}
