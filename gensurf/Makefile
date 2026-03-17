.SUFFIXES: .f .o .cpp
CCFLAGS = -I../ann_0.2/include 
FCFLAGS = -O2 -debug all -g3 -r8 -extend-source
#FCFLAG=  -O0 -tpp7 -cm -check all -W0 -axW -r8 -traceback -extend_source -FI -I /usr/local/intel/mpich/include -Vaxlib
ANNLIB = -L../ann_0.2/lib -lANN

FORLIB= -lm
# CPPLIB = -L./ann_0.2/lib -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3/../../..  -lANN -lstdc++ -lm -lgcc -lc -lgcc 

CPPLIB = -L../../ann_0.2/lib -lANN -lstdc++ -lm 

# OBJS =   gensurf.v3.0.acc.o levsets.v3.1.hh.jreblock.o outputfiles.v2.1.o walldist.o 
OBJS = soliddevider.o gensurf.v3.0.acc.o walldist.o 

# g++    $(OBJS) $(MPILIB) $(ANNLIB) $(FORLIB)    -o b.out
all: main.gensurf0.o
	ifort -o main -debug main.gensurf0.o $(OBJS) $(CPPLIB) $(FORLIB) $(ANNLIB)

main.gensurf0.o: walldist.o soliddevider.o gensurf.v3.0.acc.o main.gensurf0.f 
	ifort $(FCFLAGS) -c main.gensurf0.f

soliddevider.o: soliddevider.f90
	ifort $(FCFLAGS) -c soliddevider.f90

gensurf.v3.0.acc.o: gensurf.v3.0.acc.f
	ifort $(FCFLAGS) -c gensurf.v3.0.acc.f

walldist.o: walldist.cpp
	g++ walldist.cpp -c $(CCFLAGS) $(ANNLIB)

exe:
	./main

clean:
	rm *.o *.dat; rm Output-Tec/*; rm output/*; rm main; rm devidedsolid/*;
