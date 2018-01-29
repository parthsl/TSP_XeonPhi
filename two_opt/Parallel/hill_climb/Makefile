CC = gcc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp
OPT = -O3
objects := 2opt
DEP = src/two_opt.c src/vnn.c src/util.c

all: 2opt

debug: CFLAGS += -g -D DEBUG
debug: 2opt

2opt: lib/hill_climb.h $(DEP) src/hill_climb.c
	$(CC) -o $@ $? $(CFLAGS) -I./lib/ -Wall $(OPT)

clean:
	rm -f $(objects)

