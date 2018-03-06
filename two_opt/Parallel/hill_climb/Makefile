CC = gcc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp
OPT = -O3
objects := 2opt
DEP = src/two_opt.c src/vnn.c src/util.c
FLAGS= 
MASSLIB_PATH=/opt/ibm/xlmass/8.1.6

ifeq ($(CC),xlc)
	FLAGS+=-lmassvp8 -lmass_simdp8 -qarch=pwr8
else 
ifeq ($(CC),icc)
	FLAGS+=-fvar-tracking -axMIC-AVX512 -mkl -fimf-precision=high -mtune=knl -xHost -qopt-prefetch=2 -funroll-loops -ffast-math 
	OPT=-Ofast
endif
endif

ifeq ($(shell uname -p),ppc64le)
	FLAGS+=-L$(MASSLIB_PATH)/lib -I$(MASSLIB_PATH)/include -lmassvp8 -lmass_simdp8 -lmass -mcpu=power8 -funroll-loops -ffast-math -mrecip=all -mveclibabi=mass
endif

all: 2opt

debug: CFLAGS += -g -D DEBUG
debug: 2opt

2opt: $(DEP) src/hill_climb.c
	$(CC) -o $@ $? $(FLAGS) -g -I./lib/ -Wall $(OPT) $(CFLAGS)

clean:
	rm -f $(objects)

