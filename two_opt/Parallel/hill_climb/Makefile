CC = gcc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp
OPT = -O3
objects := 2opt
DEP = src/two_opt.c src/vnn.c src/util.c
ifeq ($(CC),xlc)
	OPT+=-lmassvp8 -lmass_simdp8 -qarch=pwr8
endif

all: 2opt

debug: CFLAGS += -g -D DEBUG
debug: 2opt

2opt: $(DEP) src/hill_climb.c
	$(CC) -o $@ $? $(CFLAGS) -g -I./lib/ -Wall $(OPT)

clean:
	rm -f $(objects)

