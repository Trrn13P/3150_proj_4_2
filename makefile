CPPflags=  clang++ -Xpreprocessor -fopenmp -std=c++14 -Rpass=loop-vectorize
LIB = -larmadillo -llapack -lblas -lomp

CONDITIONS = Lattice 60 1000000 2.0 2.3 0.01 Ordered

all: compile execute

compile:
	${CPPflags} main.cpp IsingModel.cpp -o ./main.out ${LIB}

execute:
	./main.out ${CONDITIONS}
	#python3 Energy_print.py ${CONDITIONS}
	#python3 prob_calc.py ${CONDITIONS}
