CPPflags=  clang++ -Xpreprocessor -fopenmp -std=c++14 -Rpass=loop-vectorize
LIB = -larmadillo -llapack -lblas -lomp

CONDITIONS = Lattice 2 1000000 1.0 1.1 0.1 Ordered

all: compile execute

compile:
	${CPPflags} main.cpp IsingModel.cpp -o ./main.out ${LIB}

execute:
	./main.out ${CONDITIONS}
	#python3 Energy_print.py ${CONDITIONS}
