include example.mk

CC=mpic++

CFLAGS= -O0
#CFLAGS= -O2
LFLAGS= -O0
#LFLAGS= -O2

LDIR =
RUNDIR = ../wrk/

OBJ = main.o kernel.o calculations.o test.o

#libary rules
$(OBJ): global.h findNeighbors.h 

%.o: %.cpp
	$(CC) $(CFLAGS) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

stoic: $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS_PATH) $(LIBS)
	mv stoic $(RUNDIR)

all: stoic

prepare: 
	mkdir -p $(RUNDIR)

run: all
	mpirun --oversubscribe -np 2 $(RUNDIR)/stoic

debug:
	../../../gdbgui/launch_mpi_debugger 2 $(RUNDIR)/stoic

.PHONY: clean all run

clean:
	rm -f *.o *~ core $(RUNDIR)stoic
	rm -f $(RUNDIR)*.vtp $(RUNDIR)*.pvtp

flags:
	echo "RUNDIR=" $(RUNDIR) 
	echo "CFLAGS=" $(CFLAGS) 
	echo "LFLAGS=" $(LFLAGS) 

