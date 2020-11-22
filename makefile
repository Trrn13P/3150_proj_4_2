CPPflags=  clang++ -Xpreprocessor -fopenmp -std=c++14 -Rpass=loop-vectorize
LIB = -larmadillo -llapack -lblas -lomp

CONDITIONS = Lattice 5 100000 2.0 2.3 0.1 Ordered

all: compile execute

compile:
	${CPPflags} main.cpp IsingModel.cpp -o ./main.out ${LIB}

execute:
	./main.out ${CONDITIONS}
	#python3 ./python_files/Energy_print.py ${CONDITIONS}
	#python3 ./python_files/prob_calc.py ${CONDITIONS}
	#python3 ./python_files/paralleization_analysation
	#python3 ./python_files/lattice_reader.py
