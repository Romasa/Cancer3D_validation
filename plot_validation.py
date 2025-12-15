import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#Plotting points from the experiment data\
    
    #EGFR Overexpression
xover = [7.040287769784173, 14.010071942446043 ]
yover = [68.6084142394822, 160.5177993527508]

    #EGFR Normal
x_dmso = [7.040287769784173, 14.010071942446043]
y_dmso = [49.190938511326856, 89.32038834951456]


#Reading and plotting points from simulation data
file = "EGF0.4"
file_type = ".csv"
diff = 24*60/0.02
day_ticks = np.arange(0*diff, 15*diff, diff)
xa = np.arange(1, 16, 1)

colnames = ['time', 'cell_count', 'egf_count', 'drug_count','all_erk', 'egfr_active', 'egfr_inactive', 'egfr_blocked']

folder = "06receptors/base"
dfn = pd.read_csv(folder+ '/' + "data.txt", header=None, names=colnames, sep=' ')
dfn.columns = colnames
dfn_days = dfn[dfn['time'].isin(day_ticks)]
yn = dfn_days['cell_count']
yn = yn.to_list()
yn.insert(0, 1)

folder = "24receptors/base"
dfo = pd.read_csv(folder+ '/' + "data.txt", header=None, names = colnames, sep=' ')
dfo_days = dfo[dfo['time'].isin(day_ticks)]

yh_over = dfo_days['cell_count']
yh_over = yh_over.to_list()
yh_over.insert(0,1)

def read_data(filename):
    with open(filename, mode='r', newline='') as file:
        reader = csv.reader(file)
        data_list = list(reader)
    
    return data_list
    
y_cont = read_data('cell_count_6_receptors.csv') #Normal with 6 receptor clusters
y_cont_n = [float(item) for item in y_cont[0]]

y_cont_o = read_data('cell_count_24_receptors.csv') # Overexpression with 24 receptor clusters
y_cont_over = [float(item) for item in y_cont_o[0]]

#Plotting the figures

fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.set_xlabel('Time (days)', fontsize=14)


# Plotting the result of hybrid model
ax1.set_ylabel('Simulated cell numbers', color='teal', fontsize=14)
ax1.tick_params(axis='y', colors='teal', labelsize=14)

ax1.plot(xa, yn, label="Hybrid model w/o oe-EGFR", color='teal', linewidth=2, linestyle='solid')

ax1.plot(xa, yh_over, label="Hybrid model + oe-EGFR", color='teal', linewidth=2, linestyle='solid')

# Plotting the result of continuous model
ax1.plot(xa, y_cont_n, label='Cont. Model w/o oe-EGFR', color='purple', linewidth=2, linestyle='dotted')
ax1.plot(xa, y_cont_over, label='Cont. Model + oe-EGFR', color='purple', linewidth=2, linestyle='dotted')
ax1.legend(loc='upper left')
ax1.set_ylim([0, 125])

ax1.grid(True)
# Plotting lab-experiment results
ax2 = ax1.twinx()
ax2.set_ylabel('Tumor volume in lab experiment ($mm^3$)', color='red', fontsize=14)
ax2.tick_params(axis='y', colors='red', labelsize=14)
ax2.scatter(x_dmso, y_dmso, label='lab experiment w/o oe-EGFR', color='red', linewidth=6, marker="o")
ax2.scatter(xover, yover, label='lab experiment + oe-EGFR', color='red', linewidths=6, marker=",")
ax2.set_ylim([25, 200])
ax2.legend(loc='upper right')
ax2.grid(True)

plt.grid(which='major', axis='both')

plt.savefig('validation_combined.png', dpi=300)
plt.show()
