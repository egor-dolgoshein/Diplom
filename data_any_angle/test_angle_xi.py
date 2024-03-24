import matplotlib.pyplot as plt
import numpy as np

IMXi_interaction_30 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.2PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_30 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.2PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_30 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.2PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_30 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.2PI_hx0.0500_ht0.0100.txt", 'r')


IMXi_interaction_60 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.3PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_60 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.3PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_60 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.3PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_60 = open("D:\\Documents\\Egor\\study\\Diplom\\data_any_angle\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.3PI_hx0.0500_ht0.0100.txt", 'r')



data_files_30 = [IMXi_interaction_30, REXi_interaction_30, IMXi_without_interaction_30, REXi_without_interaction_30]
IMXi_interaction_data_30 = []
REXi_interaction_data_30 = []
IMXi_without_interaction_data_30 = []
REXi_without_interaction_data_30 = []

data_files_60 = [IMXi_interaction_60, REXi_interaction_60, IMXi_without_interaction_60, REXi_without_interaction_60]
IMXi_interaction_data_60 = []
REXi_interaction_data_60 = []
IMXi_without_interaction_data_60 = []
REXi_without_interaction_data_60 = []


data_values_30 = [IMXi_interaction_data_30, REXi_interaction_data_30, IMXi_without_interaction_data_30, REXi_without_interaction_data_30]
data_values_60 = [IMXi_interaction_data_60, REXi_interaction_data_60, IMXi_without_interaction_data_60, REXi_without_interaction_data_60]
for i in range(4):
    plt.figure()
    abscissa = []
    lines_30 = data_files_30[i].readlines()
    lines_60 = data_files_60[i].readlines()
    for line in lines_30:
        abscissa.append(np.log10(float(line.split(" ")[0]))) 
        data_values_30[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_60:
        data_values_60[i].append(float(line.split(" ")[1][:-1]))
    data_files_30[i].close()
    data_files_60[i].close()
    if(i == 0):
        plt.ylabel('Im(X) (interaction)')
    if(i == 1):
        plt.ylabel('Re(X) (interaction)')
    if(i == 2):
        plt.ylabel('Im(X) (noninteraction)')
    if(i == 3):
        plt.ylabel('Re(X) (noninteraction)')
    #plt.axes(xlim=(round(abscissa[0],3),round(abscissa[len(abscissa) - 1],3)), ylim=(min(data_values[i]) - 0.1, max(data_values[i]) + 0.1))
    if(i == 0 or i == 2):
        ticks = []
        for j in range(13):ticks.append(j/20.)
    else:
        ticks = [-0.02]
        for j in range(25):ticks.append(j/20.) #[0, 1,2]
    plt.yticks(ticks)
    plt.xlabel('wt')
    plt.plot(abscissa, data_values_30[i], color = 'blue', marker="o", markersize=3, linestyle = 'dashed', label='Numerical approximation with angle 30 deg')
    plt.plot(abscissa, data_values_60[i], color = 'green', marker="o", markersize=3, linestyle = 'dashed', label='Numerical approximation with angle 60 deg')
    #plt.plot(abscissa, X_val, color='green', marker='o', markersize=3., label='Exact solution')
    plt.legend()
    plt.show()

