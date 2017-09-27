CC = gcc
CFLAGS = -fopenmp -lm
objects := two_opt

two_opt : two_opt.c
	$(CC) -o $@ $? $(CFLAGS)
clean:
	rm *.o $(objects)
