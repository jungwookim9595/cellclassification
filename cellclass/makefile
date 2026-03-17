FC = mpif90
TARGET = cellclass.exe

FCFLAGS = -O2 -debug all -g3 -r8 -extend-source
CCFLAGS = -I../ann_0.2/include

ANNLIB = -L../ann_0.2/lib -lANN
FORLIB= -lm
CPPLIB = -L../ann_0.2/lib -lANN -lstdc++ -lm 

global      = module_global.o
module_mpi  = module_mpi_topology.o module_mpi_subdomain.o
walldist	= walldist.o
cellclass 	= module_mpi_cellclass.o
main 	    = main.o 

OBJS   	    = $(global) $(module_mpi) $(walldist) $(cellclass) $(main)

all : main.o
	$(FC) -o $(TARGET) $(OBJS) $(CPPLIB) $(FORLIB) $(ANNLIB)

module_global.o: module_global.f90
	$(FC) -c module_global.f90

module_mpi_topology.o: module_mpi_topology.f90
	$(FC) -c module_mpi_topology.f90

module_mpi_subdomain.o: module_mpi_subdomain.f90
	$(FC) -c module_mpi_subdomain.f90

walldist.o: walldist.cpp
	mpicxx walldist.cpp -c $(CCFLAGS) $(ANNLIB)

module_mpi_cellclass.o: walldist.o module_mpi_cellclass.f90
	$(FC) -c module_mpi_cellclass.f90

main.o: $(global) $(module_mpi) $(walldist) $(cellclass) main.f90
	$(FC) -c main.f90

clean : 
	rm a.out *.mod *.o output/debug/* output/result/G/* output/result/SDF/*

exe :
	mpirun -np 16 ./cellclass.exe >> output/runnung_log.txt

