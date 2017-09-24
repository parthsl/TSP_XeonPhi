#include <omp.h>
#include <stdio.h>

#define numthreads 4


const int strips = 100000;

int main(int argc, int* argv){
	double st = omp_get_wtime();
	double sum=0, stripwidth=1.0/(double)strips;

	#pragma omp parallel
	#pragma omp critical
	for(int i=0;i<strips;i++){
		double x = (double)i*stripwidth;
		sum = sum + (4/(1+x*x));
	}

	printf("Pi = %lf\n",sum/strips);

	printf("Time Taken = %g\n",omp_get_wtime()-st);
	return 0;
}
