# ==============================================================================
# GenSurf Makefile - Updated for nvfortran build environment
# Converted from ifort to nvfortran (2026-03-23)
# ==============================================================================

.SUFFIXES: .f90 .o .cpp

# Compilers
FC = nvfortran
CXX = g++

# Fortran flags
# -r8 flag retained as backup for any remaining default-real code paths.
# Primary precision control is via explicit real(8) in source code.
# Remove -r8 once all real declarations are confirmed as real(8).
# FCFLAGS = -O2
# Debug build:
FCFLAGS = -O0 -g -Mbounds -Mchkptr -Minfo

# C++ flags
CCFLAGS = -I../ann_0.2/include

# Libraries
ANNLIB = -L../ann_0.2/lib -lANN
CPPLIB = -L../ann_0.2/lib -lANN -lstdc++ -lm
FORLIB = -lm

# Object files
OBJS = ann_interface.o gensurf.o soliddevider.o walldist.o

# ==============================================================================
# Targets
# ==============================================================================

all: main

main: main_gensurf.o $(OBJS)
	$(FC) $(FCFLAGS) -o main main_gensurf.o $(OBJS) $(CPPLIB) $(FORLIB) $(ANNLIB)

# Module dependencies: gensurf.o provides the GenSurf module
# ann_interface.o provides the ann_interface module
# main_gensurf.o depends on GenSurf module
# gensurf.o depends on GenSurf module (self-contained)

main_gensurf.o: main_gensurf.f90 gensurf.o
	$(FC) $(FCFLAGS) -c main_gensurf.f90

gensurf.o: gensurf.f90
	$(FC) $(FCFLAGS) -c gensurf.f90

ann_interface.o: ann_interface.f90
	$(FC) $(FCFLAGS) -c ann_interface.f90

soliddevider.o: soliddevider.f90
	$(FC) $(FCFLAGS) -c soliddevider.f90

walldist.o: walldist.cpp
	$(CXX) $(CCFLAGS) -c walldist.cpp

# ==============================================================================
# Utility targets
# ==============================================================================

exe:
	./main

clean:
	rm -f *.o *.mod *.dat main
	rm -f Output-Tec/* output/* devidedsolid/*

.PHONY: all exe clean
