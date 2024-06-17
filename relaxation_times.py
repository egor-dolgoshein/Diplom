import matplotlib.pyplot as plt
import numpy as np

data_file_relax_times_sigma_1 = open("D:\\Documents\\Egor\\study\\Diplom\\Relaxation_times_sigma_1.txt", 'r')
data_file_relax_times_sigma_3 = open("D:\\Documents\\Egor\\study\\Diplom\\Relaxation_times_sigma_3.txt", 'r')

relax_times_1 = []
relax_times_3 = []
xi = []

for line in data_file_relax_times_sigma_1:
    xi.append(line.split()[0])
    relax_times_1.append(1./float(line.split()[1]))
for line in data_file_relax_times_sigma_3:
    relax_times_3.append(1./float(line.split()[1]))
plt.figure()
plt.ylabel("Relaxation times")
plt.xlabel("xi")
plt.xticks([0, 30, 45, 60, 90])
plt.plot(xi, relax_times_1,color="blue", marker='o',markersize=2, label="sigma = 1")
plt.plot(xi, relax_times_3,color="green", marker='o',markersize=2, label="sigma = 3")
plt.legend()
plt.show()