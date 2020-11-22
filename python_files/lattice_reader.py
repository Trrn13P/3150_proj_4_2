import numpy as np
import matplotlib.pyplot as plt

"""
The file format here is the same for every line
Temperature Energy Heat_capacity Suceptibilty abs(M)
"""

#Critical temperature from Lars Onsager
T_C = 2.269

#Reading from infile and returning an array for temperature, E, C, chi and abs(M).
def infile_reader(filename):
    infile = open("../textfiles/"+filename,"r")
    T_ = []
    E_ = []
    C_ = []
    chi_ = []
    Mabs_ = []
    for line in infile:
        line = line.split()
        T_.append(eval(line[0]))
        E_.append(eval(line[1]))
        C_.append(eval(line[2]))
        chi_.append(eval(line[4]))
        Mabs_.append(eval(line[5]))
    infile.close()
    return T_, E_, C_, chi_, Mabs_


filenames = ["Lattice40","Lattice60","Lattice80","Lattice100"]
num_files = len(filenames)

#Setting up arrays for holding the values for each filename.
T = [[],[],[],[]]
E = [[],[],[],[]]
C = [[],[],[],[]]
chi = [[],[],[],[]]
Mabs = [[],[],[],[]]
#Going trough each filename and saving to arrays
i = 0
for filename in filenames:
    T[i], E[i], C[i], chi[i], Mabs[i] = infile_reader(filename)
    i+=1

#Plot function with type as argument, type can be "E", "C", "chi" or "Mabs"
#This will plot these values to file repsecitivly.
def plotter(type):
    plt.axvline(x=T_C,linestyle="--",label=r"Theoretical $T_C=2.269$",color="black")
    plt.xlabel("T",FontSize=12)
    if type=="E":
        plt.ylabel(r"$\langle E\rangle$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],E[i],label=name)
        plt.legend()
        plt.savefig("../figs/f/E.png")
        plt.clf()

    elif type=="C":
        plt.ylabel(r"$C_V/T^2$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],C[i],label=name)
        plt.legend()
        plt.savefig("../figs/f/C_V.png")
        plt.clf()


    elif type=="chi":
        plt.ylabel(r"$\chi$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],chi[i],label=name)
        plt.legend()
        plt.savefig("../figs/f/chi.png")
        plt.clf()

    elif type=="Mabs":
        plt.ylabel(r"$\langle |M|\rangle$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],Mabs[i],label=name)
        plt.legend()
        plt.savefig("../figs/f/Mabs.png")
        plt.clf()

#Setting up plot types and saving plots
plot_types = ["E","C","chi","Mabs"]
for plot_ in plot_types:
    plotter(plot_)

#printing the critical temperature for each of the lattices to terminal
#The values are provided in the results of the report.
T_vals = []
for C_ in np.asarray(C):
    T_vals.append(T[-1][np.argmax(C_)])
print(T_vals)
