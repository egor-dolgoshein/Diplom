import matplotlib.pyplot as plt
import numpy as np
#xi = 90

IMXi_interaction_per_psi_0 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_per_psi_0 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_per_psi_0 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_per_psi_0 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')

IMXi_interaction_per_psi_90 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.500PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_per_psi_90 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.500PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_per_psi_90 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.500PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_per_psi_90 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.500PI_hx0.0500_ht0.0100.txt", 'r')

data_files_per_psi_0 = [IMXi_without_interaction_per_psi_0, REXi_without_interaction_per_psi_0, IMXi_interaction_per_psi_0, REXi_interaction_per_psi_0]
data_files_per_psi_90 = [IMXi_without_interaction_per_psi_90, REXi_without_interaction_per_psi_90, IMXi_interaction_per_psi_90, REXi_interaction_per_psi_90]

IMXi_interaction_data_per_psi_0 = []
REXi_interaction_data_per_psi_0 = []
IMXi_without_interaction_data_per_psi_0 = []
REXi_without_interaction_data_per_psi_0 = []

IMXi_interaction_data_per_psi_90 = []
REXi_interaction_data_per_psi_90 = []
IMXi_without_interaction_data_per_psi_90 = []
REXi_without_interaction_data_per_psi_90 = []

data_values_per_psi_0 = [IMXi_without_interaction_data_per_psi_0, REXi_without_interaction_data_per_psi_0, IMXi_interaction_data_per_psi_0, REXi_interaction_data_per_psi_0]
data_values_per_psi_90 = [IMXi_without_interaction_data_per_psi_90, REXi_without_interaction_data_per_psi_90, IMXi_interaction_data_per_psi_90, REXi_interaction_data_per_psi_90]

for i in range(4):
    plt.figure()
    abscissa = []
    lines_per_psi_0 = data_files_per_psi_0[i].readlines()
    lines_per_psi_90 = data_files_per_psi_90[i].readlines()

    for line in lines_per_psi_0:
        abscissa.append(np.log10(float(line.split(" ")[0]))) 
        data_values_per_psi_0[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_per_psi_90:
        data_values_per_psi_90[i].append(float(line.split(" ")[1][:-1]))
    
    data_files_per_psi_0[i].close()
    data_files_per_psi_90[i].close()
    
    if(i == 0):
        plt.ylabel('Im(X) (noninteraction)')
    if(i == 1):
        plt.ylabel('Re(X) (noninteraction)')
    if(i == 2):
        plt.ylabel('Im(X) (interaction)')
    if(i == 3):
        plt.ylabel('Re(X) (interaction)')
    #plt.axes(xlim=(round(abscissa[0],3),round(abscissa[len(abscissa) - 1],3)), ylim=(min(data_values[i]) - 0.1, max(data_values[i]) + 0.1))
    if(i == 0 or i == 2):
        ticks = []
        for j in range(13):ticks.append(j/20.)
    else:
        ticks = [0]
        for j in range(25):ticks.append(j/20.) #[0, 1,2]
    #plt.yticks(ticks)
    plt.xlabel('lg(wt)')
    plt.plot(abscissa, data_values_per_psi_0[i], color = 'blue', marker="o", markersize=5, label='Psi = 0')
    plt.plot(abscissa, data_values_per_psi_90[i], color = 'red', marker="o", markersize=3, linestyle = 'dashed', label='Psi = 90')
    
    #plt.plot(abscissa, X_val, color='green', marker='o', markersize=3., label='Exact solution')
   # for j in range(len(data_values_par[i])):
   #     print('original/mine = ', data_values_par[i][j] / data_values_zero_deg[i][j])
  #  print('\n')
    plt.legend()
    plt.show()
