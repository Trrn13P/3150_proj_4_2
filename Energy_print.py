import numpy as np
import matplotlib.pyplot as plt
import sys

"""
The energy-files from c++ have the following format:
First line-> C_cycles=10000 T=2.4 number_of_spins=400 order=Unordered
Second line-> Energy-vector: MagnetcMoment-vector: Number of flips:
Third to last-> Energy-value MagneticMoment-value Flip-value
"""

#Getting filename from the makefile or terminal
try:
    filename = str(sys.argv[1]) + str(sys.argv[2])
except:
    print("Error")
    sys.exit(1)
infile = open("./textfiles/"+filename,"r")

first_line = infile.readline().split()
infile.readline()

#getting the values from the first line
MCcycles = eval(first_line[0].split("=")[1])
T = eval(first_line[1].split("=")[1])
N = eval(first_line[2].split("=")[1])
order = first_line[3].split("=")[1]

#setting up vectors to hold the data from file
mc_cycles_vec = np.arange(MCcycles+1)
energy_vec = np.zeros(MCcycles+1)
magnetic_moment_vec = np.zeros(MCcycles+1)
number_of_flips_vec = np.zeros(MCcycles+1)

#going trough the infile and saving values to vectors
for i in range(0,MCcycles+1):
    line = infile.readline().split()
    energy_vec[i] = eval(line[0])
    magnetic_moment_vec[i] = eval(line[1])
    number_of_flips_vec[i] = eval(line[2])
infile.close()

#function for plotting energy vs mc-cycles
def energy():
    plt.xlabel("MC cycles",FontSize=12)
    plt.ylabel(r"$\langle E \rangle$",FontSize=12)
    plt.plot(mc_cycles_vec,energy_vec/N)
    type = "energy"
    return type

#function for plotting magnetization vs mc-cycles
def magnetization():
    plt.xlabel("MC cycles",FontSize=12)
    plt.ylabel(r"$\langle M \rangle$",FontSize=12)
    plt.plot(mc_cycles_vec,magnetic_moment_vec/N)
    type = "magnetization"
    return type

#function for plotting # of flips vs mc-cycles
def flips():
    plt.xlabel("MC cycles",FontSize=12)
    plt.ylabel("# of flips",FontSize=12)
    plt.plot(mc_cycles_vec,number_of_flips_vec)
    type = "flips"
    return type

#function for saving the plot
def save_plot(type):
    pic_filename = "./figs/d/T=" + str(T) + "_" + str(order) + "_" + str(type) + "_L" + str(int(np.sqrt(N))) + ".png"
    plt.savefig(pic_filename)
    plt.clf()

#saving all plots for a spesific T-value
type = magnetization()
save_plot(type)
type = energy()
save_plot(type)
type = flips()
save_plot(type)
