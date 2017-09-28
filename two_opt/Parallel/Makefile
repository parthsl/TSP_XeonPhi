CC = gcc
CFLAGS = -fopenmp -lm
objects := omp_two_opt

omp_two_opt : omp_two_opt.c
	$(CC) -o $@ $? $(CFLAGS)
clean:
	rm *.o $(objects)
