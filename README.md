# FYS3150 Project 4

To run the program, go into the directory of the makefile and write "make all" in the terminal.

The project is mainly run in the makefile. There is a line called CONDITIONS which you can change to make the program do different things. The first argument can be "Lattice", "Energy" or "parallelization_test". If Lattice it will save the expectation values for a given lattice to file called LatticeXX. If Energy it will save the energy, magnetization and # of flops to file called EnergyXX. And If parallelization_test it will save the time used for 10 runs of N MCcycles parallelized and unparallelized to a file called parallelization_testXX.

The second argument the size of the lattice L, and what I in the examples over have put as XX. The third argument is the number of MCcycles. The fourth to sixth argument is respectivly the Initial temperature, Final temperature and the temperature step size. 

The last argument is "Ordered" or "Unordered". If Ordered all the Spins in the SpinMatrix will be set to 1.0. If Unordered all will be assigned a random value that is -1.0 or 1.0.

If you remove the comments with the python functions it also run those.

All python files will save plots, they are ment to be run from the makefile or the terminal, now they are set up to run from the terminal so you will get an error when running the python files from the makefile due to the directory of the textfiles/figures. If you want to run it from makefile you need to go into the python files and change the directory paths from "../" to "./".

The python programs:

Energy_print.py saves plots of the energy, magnetization and # of flips for a given L. To run it you need to give it two arguments that will be the filename it opens. This is done so that you can run it with the same CONDITIONS in the makefile and give it for example Energy 20 to open and plot from that file. 

prob_calc.py saves a plot of the normal distribution with a given variance from a lattice over a plot of the probability of finding the system at a given energy. This inputs the same file as Energy_print.py and therefore requires the same system arguments.

lattice_reader.py saves plots of the abs(M), Energy, suceptibility and heat capacity as a function of Temperature for all lattices. If you want to add more lattices you need to change the list called "filenames" at line 31 (after getting textfiles from running makefile). This also prints the critical temperatures for each of the lattices to terminal(this is the temperature at the max value of the heat capacity).

parallelization_analysation.py plots the mean of the 10 runs in parallelization_test for all lattices in the list "filenames" at line 14. To run more of these, add them to the filenames (after running the makefile).
