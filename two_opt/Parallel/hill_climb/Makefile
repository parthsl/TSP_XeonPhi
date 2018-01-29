CC = gcc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp -g
WALL = 
objects := 2opt
DEP = src/two_opt.c src/vnn.c src/util.c
OPT = -O3 -mcpu=power8 -mtune=power8 -funroll-loops -ffast-math
VEC = -maltivec -ftree-vectorizer-verbose=5

all: 2opt

debug: CFLAGS += -g -D DEBUG
debug: 2opt

2opt: $(DEP) src/hill_climb.c
	$(CC) -o $@ $? $(CFLAGS) -I./lib/ $(WALL) $(OPT) $(VEC)

2opt_xlc: $(DEP) src/hill_climb.c
	xlc -o $@ $? $(CFLAGS) -I./lib $(WALL)

clean:
	rm -f $(objects)

