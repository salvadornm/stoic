include example.mk

CC=mpic++
CFLAGS+= -v

LDIR =
RUNDIR = ../wrk/

OBJ = main.o kernel.o findNeighbors.o

%.o: %.cpp
	$(CC) -O0 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

stoic: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
	mv stoic $(RUNDIR)

all: stoic

prepare: 
	mkdir -p $(RUNDIR)

run: all
	mpirun --oversubscribe -np 2 ./stoic

debug:
	../../../gdbgui/launch_mpi_debugger 2 ./stoic

.PHONY: clean all run

clean:
	rm -f *.o *~ core stoic
	rm *.vtp *.pvtp

