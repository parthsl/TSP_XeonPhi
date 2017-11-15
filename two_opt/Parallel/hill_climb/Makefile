CC = gcc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp
WALL = -WALL
objects := 2opt 2opt_Pr
DEP = src/two_opt.c src/vnn.c src/util.c src/bruteforcer.c

all: 2opt

debug: CFLAGS += -g -D DEBUG
debug: 2opt

2opt: lib/hill_climb.h $(DEP) src/hill_climb.c
	$(CC) -o $@ $? $(CFLAGS) -I./lib/ -Wall

2opt_Pr: lib/hill_climb.h $(DEP) src/hill_climb_with_Pr.c
	$(CC) -o $@ $? $(CFLAGS) -I./lib/ -Wall

clean:
	rm -f $(objects)

