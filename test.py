import matplotlib.pyplot as plt
import numpy as np

#interaction off
def X_(wt):
    return 1./(1+wt*wt)

def X__(wt):
    return X_(wt)*wt

#interaction on

def X_inter(wt): return X_(wt) + (X_(wt)*X_(wt) - X__(wt)*X__(wt))/3.

def X__inter(wt): return X__(wt)*(1+2*X_(wt)/3.)

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
    X_val = []
    if(i == 0):
        plt.ylabel('Im(X) (interaction)')
        for num in abscissa: X_val.append(X__inter(np.exp(num)))
    if(i == 1):
        plt.ylabel('Re(X) (interaction)')
        for num in abscissa: X_val.append(X_inter(np.exp(num)))
    if(i == 2):
        plt.ylabel('Im(X) (noninteraction)')
        for num in abscissa: X_val.append(X__(np.exp(num)))
    if(i == 3):
        plt.ylabel('Re(X) (noninteraction)')
        for num in abscissa: X_val.append(X_(np.exp(num)))
    #plt.axes(xlim=(round(abscissa[0],3),round(abscissa[len(abscissa) - 1],3)), ylim=(min(data_values[i]) - 0.1, max(data_values[i]) + 0.1))
    if(i == 0 or i == 2):
        ticks = []
        for j in range(13):ticks.append(j/20.)
    else:
        ticks = [-0.02]
        for j in range(25):ticks.append(j/20.) #[0, 1,2]
    plt.yticks(ticks)
    plt.xlabel('wt')
    plt.plot(abscissa, data_values[i], color = 'blue', marker="o", markersize=3, linestyle = 'dashed', label='Numerical approximation')
    #plt.plot(abscissa, X_val, color='green', marker='o', markersize=3., label='Exact solution')
    plt.legend()
    plt.show()

