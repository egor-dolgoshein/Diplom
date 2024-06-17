import matplotlib.pyplot as plt
import numpy as np

IMXi_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\IMXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\IMXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\REXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_par = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\REXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')

IMXi_interaction_par_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_par_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_par_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_par_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\0_deg_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.000PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')


IMXi_interaction_xi_30 = open("D:\\Documents\\Egor\\study\\Diplom\\30_deg_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.167PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_xi_30 = open("D:\\Documents\\Egor\\study\\Diplom\\30_deg_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.167PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_xi_30 = open("D:\\Documents\\Egor\\study\\Diplom\\30_deg_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.167PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_xi_30 = open("D:\\Documents\\Egor\\study\\Diplom\\30_deg_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.167PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')


IMXi_interaction_xi_45 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\IMXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_xi_45 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\IMXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_xi_45 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\REXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_xi_45 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\REXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')

IMXi_interaction_xi_45_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_xi_45_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_xi_45_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_xi_45_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\45_deg_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.250PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')


IMXi_interaction_xi_60 = open("D:\\Documents\\Egor\\study\\Diplom\\60_deg_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.333PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_xi_60 = open("D:\\Documents\\Egor\\study\\Diplom\\60_deg_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.333PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_xi_60 = open("D:\\Documents\\Egor\\study\\Diplom\\60_deg_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.333PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_xi_60 = open("D:\\Documents\\Egor\\study\\Diplom\\60_deg_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.333PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')


IMXi_interaction_per = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_per = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_per = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_per = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_without_interaction_alpha0.0_sigma3.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')

IMXi_interaction_per_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
IMXi_without_interaction_per_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\IMXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_interaction_per_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_interaction_alpha0.0_sigma1.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')
REXi_without_interaction_per_sigma1 = open("D:\\Documents\\Egor\\study\\Diplom\\90_deg_data\\REXi_without_interaction_alpha0.0_sigma1.0_chil1.0_xi0.500PI_psi0.000PI_hx0.0500_ht0.0100.txt", 'r')


data_files_per = [IMXi_without_interaction_per, REXi_without_interaction_per, IMXi_interaction_per, REXi_interaction_per]
data_files_par = [IMXi_without_interaction_par, REXi_without_interaction_par, IMXi_interaction_par, REXi_interaction_par]
IMXi_interaction_data_par = []
REXi_interaction_data_par = []
IMXi_without_interaction_data_par = []
REXi_without_interaction_data_par = []

IMXi_interaction_data_per = []
REXi_interaction_data_per = []
IMXi_without_interaction_data_per = []
REXi_without_interaction_data_per = []


data_files_xi_45 = [IMXi_without_interaction_xi_45, REXi_without_interaction_xi_45, IMXi_interaction_xi_45, REXi_interaction_xi_45]
IMXi_interaction_data_xi_45 = []
REXi_interaction_data_xi_45 = []
IMXi_without_interaction_data_xi_45 = []
REXi_without_interaction_data_xi_45 = []

data_files_per_sigma1 = [IMXi_without_interaction_per_sigma1, REXi_without_interaction_per_sigma1, IMXi_interaction_per_sigma1, REXi_interaction_per_sigma1]
data_files_par_sigma1 = [IMXi_without_interaction_par_sigma1, REXi_without_interaction_par_sigma1, IMXi_interaction_par_sigma1, REXi_interaction_par_sigma1]
IMXi_interaction_data_par_sigma1 = []
REXi_interaction_data_par_sigma1 = []
IMXi_without_interaction_data_par_sigma1 = []
REXi_without_interaction_data_par_sigma1 = []

IMXi_interaction_data_per_sigma1 = []
REXi_interaction_data_per_sigma1 = []
IMXi_without_interaction_data_per_sigma1 = []
REXi_without_interaction_data_per_sigma1 = []


data_files_xi_45_sigma1 = [IMXi_without_interaction_xi_45_sigma1, REXi_without_interaction_xi_45_sigma1, IMXi_interaction_xi_45_sigma1, REXi_interaction_xi_45_sigma1]
IMXi_interaction_data_xi_45_sigma1 = []
REXi_interaction_data_xi_45_sigma1 = []
IMXi_without_interaction_data_xi_45_sigma1 = []
REXi_without_interaction_data_xi_45_sigma1 = []

data_files_xi_30 = [IMXi_without_interaction_xi_30, REXi_without_interaction_xi_30, IMXi_interaction_xi_30, REXi_interaction_xi_30]
IMXi_interaction_data_xi_30 = []
REXi_interaction_data_xi_30 = []
IMXi_without_interaction_data_xi_30 = []
REXi_without_interaction_data_xi_30 = []

data_files_xi_60 = [IMXi_without_interaction_xi_60, REXi_without_interaction_xi_60, IMXi_interaction_xi_60, REXi_interaction_xi_60]
IMXi_interaction_data_xi_60 = []
REXi_interaction_data_xi_60 = []
IMXi_without_interaction_data_xi_60 = []
REXi_without_interaction_data_xi_60 = []

data_values_per = [IMXi_without_interaction_data_per, REXi_without_interaction_data_per, IMXi_interaction_data_per, REXi_interaction_data_per]
data_values_par = [IMXi_without_interaction_data_par, REXi_without_interaction_data_par, IMXi_interaction_data_par, REXi_interaction_data_par]
data_values_xi_45 = [IMXi_without_interaction_data_xi_45, REXi_without_interaction_data_xi_45, IMXi_interaction_data_xi_45, REXi_interaction_data_xi_45]
data_values_xi_30 = [IMXi_without_interaction_data_xi_30, REXi_without_interaction_data_xi_30, IMXi_interaction_data_xi_30, REXi_interaction_data_xi_30]
data_values_xi_60 = [IMXi_without_interaction_data_xi_60, REXi_without_interaction_data_xi_60, IMXi_interaction_data_xi_60, REXi_interaction_data_xi_60]

data_values_per_sigma1 = [IMXi_without_interaction_data_per_sigma1, REXi_without_interaction_data_per_sigma1, IMXi_interaction_data_per_sigma1, REXi_interaction_data_per_sigma1]
data_values_par_sigma1 = [IMXi_without_interaction_data_par_sigma1, REXi_without_interaction_data_par_sigma1, IMXi_interaction_data_par_sigma1, REXi_interaction_data_par_sigma1]
data_values_xi_45_sigma1 = [IMXi_without_interaction_data_xi_45_sigma1, REXi_without_interaction_data_xi_45_sigma1, IMXi_interaction_data_xi_45_sigma1, REXi_interaction_data_xi_45_sigma1]

for i in range(4):
    plt.figure()
    abscissa = []
    lines_per = data_files_per[i].readlines()
    lines_par = data_files_par[i].readlines()
    lines_xi_45 = data_files_xi_45[i].readlines()
    lines_xi_30 = data_files_xi_30[i].readlines()
    lines_xi_60 = data_files_xi_60[i].readlines()
    lines_per_sigma1 = data_files_per_sigma1[i].readlines()
    lines_par_sigma1 = data_files_par_sigma1[i].readlines()
    lines_xi_45_sigma1 = data_files_xi_45_sigma1[i].readlines()

    for line in lines_par:
        abscissa.append(np.log10(float(line.split(" ")[0]))) 
        data_values_par[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_per:
        data_values_per[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_xi_45:
        data_values_xi_45[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_xi_30:
        data_values_xi_30[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_xi_60:
        data_values_xi_60[i].append(float(line.split(" ")[1][:-1]))
        
    for line in lines_par_sigma1:
        data_values_par_sigma1[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_per_sigma1:
        data_values_per_sigma1[i].append(float(line.split(" ")[1][:-1]))
    for line in lines_xi_45_sigma1:
        data_values_xi_45_sigma1[i].append(float(line.split(" ")[1][:-1]))

    data_files_par[i].close()
    data_files_per[i].close()
    data_files_xi_45[i].close()
    data_files_xi_30[i].close()
    data_files_xi_60[i].close()
    
    data_files_par_sigma1[i].close()
    data_files_per_sigma1[i].close()
    data_files_xi_45_sigma1[i].close()
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
    #plt.plot(abscissa, data_values_par[i], color = 'blue', marker="o", markersize=3, linestyle = 'dashed', label='Xi = 0')
    plt.plot(abscissa, data_values_xi_30[i], color = 'black', marker="o", markersize=3, linestyle = 'dashed', label='Xi = 30')
    #plt.plot(abscissa, data_values_xi_45[i], color = 'green', marker="o", markersize=3, linestyle = 'dashed', label='Xi = 45')
    plt.plot(abscissa, data_values_xi_60[i], color = 'yellow', marker="o", markersize=3, linestyle = 'dashed', label='Xi = 60')
    #plt.plot(abscissa, data_values_per[i], color = 'red', marker="o", markersize=3, linestyle = 'dashed', label='Xi = 90')
    plt.plot(abscissa, data_values_par_sigma1[i], color = 'blue', marker="o", markersize=3, label='Xi = 0')
    plt.plot(abscissa, data_values_xi_45_sigma1[i], color = 'green', marker="o", markersize=3, label='Xi = 45')
    plt.plot(abscissa, data_values_per_sigma1[i], color = 'red', marker="o", markersize=3, label='Xi = 90')
    
    #plt.plot(abscissa, X_val, color='green', marker='o', markersize=3., label='Exact solution')
   # for j in range(len(data_values_par[i])):
   #     print('original/mine = ', data_values_par[i][j] / data_values_zero_deg[i][j])
  #  print('\n')
    plt.legend()
    plt.show()

