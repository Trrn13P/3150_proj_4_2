import numpy as np
import matplotlib.pyplot as plt


T_C = 2.269

def infile_reader(filename):
    infile = open("./textfiles/"+filename,"r")
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

T = [[],[],[],[]]
E = [[],[],[],[]]
C = [[],[],[],[]]
chi = [[],[],[],[]]
Mabs = [[],[],[],[]]
i = 0
for filename in filenames:
    T[i], E[i], C[i], chi[i], Mabs[i] = infile_reader(filename)
    i+=1


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
        plt.savefig("./figs/f/E.png")
        plt.clf()

    elif type=="C":
        plt.ylabel(r"$C_V/T^2$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],C[i],label=name)
        plt.legend()
        plt.savefig("./figs/f/C_V.png")
        plt.clf()


    elif type=="chi":
        plt.ylabel(r"$\chi$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],chi[i],label=name)
        plt.legend()
        plt.savefig("./figs/f/chi.png")
        plt.clf()

    elif type=="Mabs":
        plt.ylabel(r"$\langle |M|\rangle$",FontSize=12)
        for i in range(num_files):
            lattice_number = filenames[i][7:]
            name = lattice_number + "x" + lattice_number
            plt.plot(T[0],Mabs[i],label=name)
        plt.legend()
        plt.savefig("./figs/f/Mabs.png")
        plt.clf()


plot_types = ["E","C","chi","Mabs"]
for plot_ in plot_types:
    plotter(plot_)
