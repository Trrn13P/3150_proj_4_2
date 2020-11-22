import numpy as np
import matplotlib.pyplot as plt

"""
There are 10 runs for each lattice, we will take the mean and plot it

The files have the format
first line: MC-cycles L
second line: t-unparellized t-paralleized
third to last: values values
"""

#setting up lists for the times
filenames = ["parallelization_test2","parallelization_test5","parallelization_test10","parallelization_test20"]
paralleized_times = [[],[],[],[]]
unparellized_times = [[],[],[],[]]

#reading from file and saving to lists
N = len(filenames)
for i in range(N):
    filename = "../textfiles/"+filenames[i]
    infile = open(filename, "r")
    infile.readline()
    infile.readline()
    for line in infile:
        line = line.split()
        unparellized_times[i].append(eval(line[0]))
        paralleized_times[i].append(eval(line[1]))

#these arrays are for holding the mean values
plot_para = []
plot_unpara = []

#getting the mean values for the times
for i in range(N):
    unparellized_time = unparellized_times[i]
    paralleized_time = paralleized_times[i]
    n = len(unparellized_time)
    unparellized_time = sum(unparellized_time)/n
    paralleized_time = sum(paralleized_time)/n
    plot_para.append(paralleized_time)
    plot_unpara.append(unparellized_time)

#saving plot
plot_x = [2,5,10,20]
plt.xlabel("L",FontSize=12)
plt.ylabel("Time [s]",FontSize=12)
plt.plot(plot_x,plot_para,label="Parallelized")
plt.plot(plot_x,plot_unpara,label="Unparallelized")
plt.savefig("../figs/para_test.png")
