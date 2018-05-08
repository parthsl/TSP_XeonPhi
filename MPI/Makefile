CC = mpicc
LIBDIR = lib
SRCDIR = src
CFLAGS = -lm -fopenmp
OPT = -O3
objects := 2opt
DEP = src/two_opt.c src/vnn.c src/util.c
FLAGS=
MASSLIB_PATH=/opt/ibm/xlmass/8.1.6
OpenMPI_PATH=/root/Parth/openmpi-3.0.0/build/
CUSTOM_OpenMPI_PATH?=no

ifeq ($(CC),mpicc)
ifeq ($(CUSTOM_OpenMPI_PATH),yes)
	FLAGS+=-I$(OpenMPI_PATH)/include -L$(OpenMPI_PATH)/lib
endif
endif
ifeq ($(CC),xlc)
	FLAGS+=-lmassvp8 -lmass_simdp8 -qarch=pwr8 -mcpu=power8 -funroll-loops
endif
ifeq ($(shell uname -p),ppc64le)
	FLAGS+=-L$(MASSLIB_PATH)/lib -I$(MASSLIB_PATH)/include -lmassvp8 -lmass_simdp8 -lmass -mcpu=power8 -funroll-loops -ffast-math -mrecip=all -mveclibabi=mass
endif

all: 2opt

debug: CFLAGS += -g -D DEBUG
debug: 2opt

2opt: $(DEP) src/hill_climb_mpich.c
	$(CC) -o $@ $? $(FLAGS) -g -I./lib/ -Wall -Wextra $(OPT) $(CFLAGS)

clean:
	rm -f $(objects)

