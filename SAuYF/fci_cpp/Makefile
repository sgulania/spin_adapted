CXX= g++
#PROFILER_FLAG=-pg 
WARNING_FLAG=-Wall -Wno-deprecated -Wfatal-errors
DEBUG_FLAG=-g  -gdwarf-2 
OPT_FLAG=  -O3 
CXXFLAGS= -std=c++11 

CI_SRC=/home/jdw/Desktop/MolecularEnsembles/tools/ci_matrix

#Sufficent libraries and includes needed to compile against libint2
#INCLUDESEig	=
#-IEigen/  -I. 

INCLUDESLibint	=-I/usr/local/include -I/usr/local/include/libint2 -I/opt/local/include
#-Ilibint/include/ -Ilibint/include/libint2

LIBSLibint      =-L/usr/local/lib -lint2
LIBSLapack	=-L/opt/local/lib/lapack -llapack  -lblas -llapacke 
LIBS	        =-lm -lgmp -lgmpxx $(LIBSLapack) $(LIBSLibint)
INCLUDES        =-I$(CI_SRC) $(INCLUDESLibint) 
#$(INCLUDESEig)

all: spin_half 
	

#test prim wecker_ensemble single_molecule_ensemble fci integral_analysis hf spin_half

spin_half.o: spin_half.cpp
	$(CXX) spin_half.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o spin_half.o

spin_half: spin_half.o ci_matrix.o parser.o libint_interface.o
	$(CXX) $(CXXFLAGS) spin_half.o libint_interface.o ensembles.o parser.o ci_matrix.o $(LIBS) -o spin_half 

test: libint_interface.o parser.o ci_matrix.o minimal.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBSLibint) -L/opt/local/lib/lapack -llapacke -llapack ci_matrix.o parser.o libint_interface.o minimal.o -o test 

parser.o: parser.cpp
	$(CXX) parser.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o parser.o

libint_interface.o: libint_interface.cpp
	$(CXX) libint_interface.cpp $(INCLUDES) -std=c++11 -c -o libint_interface.o

minimal.o: minimal.cpp
	$(CXX) minimal.cpp $(INCLUDES) -std=c++11 -c -o minimal.o


ci_matrix.o: ci_matrix.cpp 
	$(CXX) ci_matrix.cpp $(CXXFLAGS) $(INCLUDESLibint) -c -o ci_matrix.o

convert_g09: convert_g09.cpp parser.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS)  convert_g09.cpp parser.o -o convert_g09

wecker_ensemble:wecker_ensemble.o ci_matrix.o parser.o ensembles.o
	$(CXX) $(CXXFLAGS) wecker_ensemble.o ensembles.o ci_matrix.o parser.o $(LIBS) -o wecker_ensemble

fci:fci.o ci_matrix.o parser.o libint_interface.o
	$(CXX) $(CXXFLAGS) fci.o libint_interface.o ci_matrix.o ensembles.o parser.o $(LIBS) -o fci

single_molecule_ensemble:single_molecule_ensemble.o ci_matrix.o parser.o libint_interface.o
	$(CXX) $(CXXFLAGS) single_molecule_ensemble.o libint_interface.o ci_matrix.o ensembles.o parser.o $(LIBS) -o single_molecule_ensemble

single_molecule_ensemble.o: single_molecule_ensemble.cpp
	$(CXX) single_molecule_ensemble.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o single_molecule_ensemble.o

o_test: fci atomic_ensemble/o_example/basis_funcs atomic_ensemble/o_example/nuc_field
	./fci atomic_ensemble/o_example/basis_funcs atomic_ensemble/o_example/nuc_field > o_test

integral_analysis.o:integral_analysis.cpp
	$(CXX) integral_analysis.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o integral_analysis.o

integral_analysis: parser.o integral_analysis.o ci_matrix.o libint_interface.o ensembles.o
	$(CXX) $(CXXFLAGS) integral_analysis.o ensembles.o ci_matrix.o libint_interface.o  parser.o $(LIBS) -o integral_analysis

hf.o:hf.cpp
	$(CXX) hf.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o hf.o

hf: parser.o hf.o ci_matrix.o ensembles.o libint_interface.o
	$(CXX) $(CXXFLAGS) hf.o libint_interface.o ci_matrix.o parser.o ensembles.o $(LIBS) -o hf

ensembles.o: ensembles.cpp
	$(CXX) ensembles.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o ensembles.o

wecker_ensemble.o: wecker_ensemble.cpp
	$(CXX) wecker_ensemble.cpp $(CXXFLAGS) $(INCLUDES) -c -o wecker_ensemble.o


fci.o: fci.cpp ci_matrix.cpp parser.cpp
	$(CXX) fci.cpp $(CXXFLAGS) $(INCLUDESLibint) $(INCLUDESEig) -c -o fci.o


