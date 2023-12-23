import matplotlib.pyplot as plt
import numpy as np

IMXi_interaction = open("IMXi_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction = open("IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction = open("REXi_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction = open("REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')

data_files = [IMXi_interaction, REXi_interaction, IMXi_without_interaction, REXi_without_interaction]
IMXi_interaction_data = []
REXi_interaction_data = []
IMXi_without_interaction_data = []
REXi_without_interaction_data = []
data_values = [IMXi_interaction_data, REXi_interaction_data, IMXi_without_interaction_data, REXi_without_interaction_data]
for i in range(4):
    plt.figure()
    abscissa = []
    lines = data_files[i].readlines()
    for line in lines:
        abscissa.append(np.log(float(line.split(" ")[0])))
        data_values[i].append(float(line.split(" ")[1][:-1]))
    data_files[i].close()
    #plt.axes(xlim=(abscissa[0],abscissa[len(abscissa) - 1]), ylim=(data_values[i][0], data_values[i][len(data_values[i])-1]))
    plt.xlabel('wt')
    plt.ylabel('Function value')
    plt.plot(abscissa, data_values[i], color = 'blue', marker="o", markersize=3, label='XL = alpha = sigma = 1')
    plt.show()

