#
# Makefile for mpi-advection example.
#

LIBPFASST ?= ../..

EXE = main.exe

include $(LIBPFASST)/Makefile.defaults
FSRC = feval.f90 \
       hooks.f90 \
       main.f90

OBJ = $(addprefix build/,$(FSRC:.f90=.o))

FFLAGS  += -g -I$(LIBPFASST)/include -fbounds-check
LDFLAGS += -g -L$(LIBPFASST)/lib -lpfasst

all: $(EXE)

VPATHS = src

include $(LIBPFASST)/Makefile.rules

build/%.o: %.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< $(OUTPUT_OPTION)


build/%.o: %.f90
	@mkdir -p build
	@mkdir -p include
	$(FC) $(FFLAGS) -c $< $(OUTPUT_OPTION)


#
# dependencies
#

build/feval.o:
build/hooks.o: build/feval.o
build/main.o:  build/hooks.o build/feval.o
build/transfer.o:  build/feval.o
