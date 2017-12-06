CC = gcc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp -g
WALL = -WALL
objects := 2opt 2opt_Pr 2opt_mpich
DEP = src/two_opt.c src/vnn.c src/util.c

all: 2opt_mpich

debug: CFLAGS += -D DEBUG
debug: 2opt

2opt_mpich: lib/hill_climb.h $(DEP) src/hill_climb_mpich.c
	$(CC) -o $@ $? $(CFLAGS) -I./lib/ -Wall

clean:
	rm -f $(objects)

