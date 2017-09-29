CC = gcc
CFLAGS = -fopenmp -lm
objects := omp_two_opt omp_parallel_twoopt_vnn

omp_two_opt : omp_two_opt.c
	$(CC) -o $@ $? $(CFLAGS)

omp_parallel_twoopt_vnn : omp_parallel_twoopt_vnn.c
	$(CC) -o $@ $? $(CFLAGS)
 
clean:
	rm *.o $(objects)
