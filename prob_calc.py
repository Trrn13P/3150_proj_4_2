import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.integrate import trapz
import sys

try:
    filename = str(sys.argv[1]) + str(sys.argv[2])
except:
    print("Error")
    sys.exit(1)

infile = open("./textfiles/"+filename,"r")

first_line = infile.readline().split()
infile.readline()

MCcycles = eval(first_line[0].split("=")[1])
T = eval(first_line[1].split("=")[1])
N = eval(first_line[2].split("=")[1])
order = first_line[3].split("=")[1]


mc_cycles_vec = np.arange(MCcycles+1)
energy_vec = np.zeros(MCcycles+1)
magnetic_moment_vec = np.zeros(MCcycles+1)
number_of_flips_vec = np.zeros(MCcycles+1)

for i in range(0,MCcycles+1):
    line = infile.readline().split()
    energy_vec[i] = eval(line[0])
    magnetic_moment_vec[i] = eval(line[1])
    number_of_flips_vec[i] = eval(line[2])
infile.close()

#starting after 10% of the cycles
start_index = int((1/10)*MCcycles)

energy_count = Counter(energy_vec[start_index:-1])
energy_count = dict(energy_count)
magneticMoment_count = Counter(magnetic_moment_vec[start_index:-1])
magneticMoment_count = dict(magneticMoment_count)

length_ = len(energy_vec[start_index:-1])

energy_values = []
probability = []

for value in energy_count:
    energy_values.append(value)
    probability.append(energy_count[value])

lists = sorted(zip(*[energy_values, probability]))
new_energy, new_probability = list(zip(*lists))
new_probability = np.asarray(new_probability)

new_probability = new_probability/trapz(new_probability,new_energy)

plt.xlabel("E",FontSize=12)
plt.ylabel("P(E)",FontSize=12)
plt.plot(new_energy,new_probability,label="Simulated values")


E = 0
E_sqrd = 0
for i in energy_vec[start_index:-1]:
    E += i;
    E_sqrd += i**2
E = E/length_
E_sqrd = E_sqrd/length_
variance = E_sqrd-E**2



mean_E = np.sum(energy_vec[start_index:-1])/length_
#probability distribution function
f = 1/(np.sqrt(2*np.pi*variance))*np.exp(-0.5*(energy_vec[start_index:-1]-mean_E)**2/variance)

lists = sorted(zip(*[energy_vec[start_index:-1],f]))
energy_vec,f = list(zip(*lists))
f = np.asarray(f)
print(variance)
plt.plot(energy_vec,f, label="PDF")
plt.legend()
plt.savefig("./figs/e/T="+str(T)+".png")
#plt.show()
plt.clf()
