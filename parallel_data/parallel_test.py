import matplotlib.pyplot as plt
import numpy as np

IMXi_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_hx0.0500_ht0.0100.txt", 'r')


IMXi_interaction_zero_deg = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.0PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_zero_deg = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.0PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_zero_deg = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.0PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_zero_deg = open("D:\\Documents\\Egor\\study\\Diplom\\parallel_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.0PI_hx0.0500_ht0.0100.txt", 'r')



data_files_par = [IMXi_interaction_par, REXi_interaction_par, IMXi_without_interaction_par, REXi_without_interaction_par]
IMXi_interaction_data_par = []
REXi_interaction_data_par = []
IMXi_without_interaction_data_par = []
REXi_without_interaction_data_par = []

data_files_zero_deg = [IMXi_interaction_zero_deg, REXi_interaction_zero_deg, IMXi_without_interaction_zero_deg, REXi_without_interaction_zero_deg]
IMXi_interaction_data_zero_deg = []
REXi_interaction_data_zero_deg = []
IMXi_without_interaction_data_zero_deg = []
REXi_without_interaction_data_zero_deg = []


data_values_par = [IMXi_interaction_data_par, REXi_interaction_data_par, IMXi_without_interaction_data_par, REXi_without_interaction_data_par]
data_values_zero_deg = [IMXi_interaction_data_zero_deg, REXi_interaction_data_zero_deg, IMXi_without_interaction_data_zero_deg, REXi_without_interaction_data_zero_deg]
for i in range(4):
    plt.figure()
    abscissa = []
    lines_par = data_files_par[i].readlines()
    lines_zero_deg = data_files_zero_deg[i].readlines()
    for line in lines_par:
        abscissa.append(np.log10(float(line.split(" ")[0]))) 
        data_values_par[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_zero_deg:
        data_values_zero_deg[i].append(float(line.split(" ")[1][:-1]))
    data_files_par[i].close()
    data_files_zero_deg[i].close()
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
    plt.plot(abscissa, data_values_par[i], color = 'blue', marker="o", markersize=3, linestyle = 'dashed', label='Original parallel case')
    plt.plot(abscissa, data_values_zero_deg[i], color = 'green', marker="o", markersize=3, linestyle = 'dashed', label='Numerical approximation with angle 0 deg')
    #plt.plot(abscissa, X_val, color='green', marker='o', markersize=3., label='Exact solution')
    plt.legend()
    plt.show()

