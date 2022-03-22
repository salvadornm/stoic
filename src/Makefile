include example.mk

CC=mpic++

#default flags
CFLAGS= 
LFLAGS= 

LDIR =
RUNDIR = ../wrk/


OBJ = main.o kernel.o calculations.o

#libary rules
$(OBJ): global.h findNeighbors.h 


#options
ifeq ($(CANTERA),1)
  CANTERA_FLG="$(shell pkg-config --cflags --libs cantera)"
  CFLAGS += $(CANTERA_FLG)
  CFLAGS += -Dchem_cantera
endif

ifeq ($(DEBUG),1)
  CFLAGS += -Ddebug
  CFLAGS += -O0
else
  CFLAGS += -O2
endif




%.o: %.cpp
	$(CC) $(CFLAGS) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

stoic: $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS_PATH) $(LIBS)
	mv stoic $(RUNDIR)
	$(call echogreen," SUCCESS")

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
	$(call echoblue," RUNDIR =  ")
	@echo $(RUNDIR)
	$(call echoblue," compile flags  CFLAGS=  ")
	@echo $(CFLAGS)
	$(call echoblue," linking flags  LFLAGS=  ")
	@echo "LFLAGS = $(LFLAGS) "

#=============== color schemes 
define echored
      @tput setaf 1
      @echo $1
      @tput sgr0
endef
define echogreen
      @tput setaf 2
      @echo $1
      @tput sgr0
endef
define echoblue
      @tput setaf 4
      @echo $1
      @tput sgr0
endef
define echocyan
      @tput setaf 6
      @echo $1
      @tput sgr0
endef
define echoinfo
      @tput setaf 2
      @echo $1 "  :  "  `tput setaf 0` $2
      @tput sgr0
endef	 
