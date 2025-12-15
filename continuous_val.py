import csv
import numpy as np
import pandas as pd

#This routine converts the concentration of EGF into the number of molecules and downsize it by the factor of 100

avagadro_no = 6.022e23 # unit : per mole
dividing_factor = 100

mw_egf = 6400 # in Daltons
vol_sol = 50 # (radius of the computational domain) in micrometers,
nano = 10**-9

def write_to_file(CT, filename):
    lst_cellcount = []
    for i in range(15):
      lst_cellcount.append(round(CT[0][i*4800],0))
    #lst_cellcount.append(round(CT[0][(15*4800)-1],0))

    with open(filename,  'w', newline='') as myfile:
         wr = csv.writer(myfile)
         wr.writerow(lst_cellcount)
    return lst_cellcount
    
    
def radius2mL(r): # r in micrometer
  vol = (4*np.pi*r**3)/3 # unit: cubic micrometer
  vol/= 10**12 # unit: milliliter
  return vol


# input: concentration 'conc' in nanograms/mL, molecular weight 'mol_wt' in grams/mole, volume of solution 'v_sol' in micrometer
def conc2count(conc, mol_wt=mw_egf, v_sol=vol_sol, unit=nano):
    n_moles = (conc*unit)/mol_wt # unit: mole/mL
    n_molecules = n_moles*avagadro_no # unit: molecules/mL
    mol_in_sol = n_molecules*radius2mL(v_sol) # unit: molecules/micrometer
    return mol_in_sol/dividing_factor


df_growth = pd.DataFrame(columns=['Nr', 'd_egfr','egf_count', 'growth'])
df_growth

def integrate_vector(c, dx):
    return np.trapz(c, dx=dx)

def rate(n, n_sat, b):
        return b*n/(1 + 1*n)

def cell_count(k4, Nr_lst, e0_lst):
  CT = [None] * len(e0_lst)
  CTT = [None] * len(e0_lst)

  for n_key, Nr in enumerate(Nr_lst):
    for k, e0 in enumerate(e0_lst):
      cT = np.zeros([int(Nd/dt), ])
      eT = np.zeros([int(Nd/dt), ])
      gT = np.zeros([int(Nd/dt), ])

      e = e0
      e_new = e0


      c = np.zeros([int(Nr/dx) + 3, ])
      n = np.linspace(-dx, Nr + dx, int(Nr/dx)+3)

      c_new = np.zeros([int(Nr/dx) + 3, ])
      cTT = np.zeros([int(Nr/dx)+3, int(Nd/dt)+1])

      eT[0] = e0

      for i in range(int(Nr/dx)):
          if (n[i] <= int(1/dx) and n[i] >= int(0/dx)):
              c[i] = 1

      cTT[:, 0] = c

      cT[0] = integrate_vector(c, dx)

      gT[0] = integrate_vector(c*n, dx)/integrate_vector(c, dx)


      for t in range(1,int(Nd/dt)):
          e_new = e + dt*(k1*(e0-e) - alpha*k2*integrate_vector((Nr-n)*c, dx)*e - k3*e)


          for i in range(1, int(Nr/dx)):
              if n[i] < 0.5 * Nr:
                  production = (
                      0.5 * rate(n[2*i], n_sat, k4) * c[2*i] * np.log(N / max(integrate_vector(c, dx), 0.1)) +
                      + 0.5 * rate(n[2*i+1], n_sat, k4) * c[2*i+1] * np.log(N / max(integrate_vector(c, dx), 0.1))

                  )

                  c_new[i] = (
                      c[i]
                      + dt * production
                      - dt * k5 * c[i]
                      + dt * k6 * n[i]*(c[i+1] - c[i]) / dx  + dt*k6*c[i]
                      - dt * k2 * e * (Nr - n[i])*(c[i] -  c[i-1]) / dx + dt*k2* e *c[i]
                  )
              else:
                  c_new[i] = (
                      c[i]
                      - dt * k5 * c[i]
                      + dt * k6 * n[i]*(c[i+1] - c[i]) / dx  + dt*k6*c[i]
                      - dt * k2 * e * (Nr - n[i])*(c[i] -  c[i-1]) / dx + dt * k2 * e *c[i]
                  )

          e = e_new
          eT[t] = e
          gT[t] = integrate_vector(c*n, dx)/integrate_vector(c, dx)
          c = c_new

          cTT[:, t] = c

          cT[t] = integrate_vector(c, dx)

      CT[k] = list(cT)
      CTT[k] = cTT
      print("Cells at end=", CT[k][-1])

      df_growth.loc[-1] = [Nr, k6, e0, CT[k][-1]]
      df_growth.index = df_growth.index + 1

  return cTT, CT

Nr = 6
Nr_lst = [Nr]
Nd = 15*24   #number of hours

dx = 0.5
dt = 0.005

egf = np.arange(0.2, 1.01, 0.2)
#time unit: 1h
#space unit: number of receptors
k1 = 50 #egf influx rate
k2 = 0.0025 #egf-egfr binding rate
k3 = 0.69  #egf degradation rate  ~ 4 hours
k4 = 1/(24*2.875) # division rate when once receptor is active, for 6 receptors
#k4 = 1/(24*3.5) # division rate when once receptor is active, for 24 receptors
k5 = 1/24 #apoptosis rate   ~ 24 hours

n_sat = 4    #saturation receptor
alpha = 1   #consumption rate
N = 5e5 # cell capacity

k6 = 1/1.5  #egfr inactivation
e0_conc = [0.05, 0.1, 0.15, 0.2, 0.4, 0.6]
e0_conc = [0.4]
e0_lst =  [round(conc2count(x),0) for x in e0_conc]    #initial egf


#postprocessing vectors:
cT = np.zeros([int(Nd/dt), ])
eT = np.zeros([int(Nd/dt), ])

cTT, CT = cell_count(1/(24*2.875) , [6], e0_lst)
cell_growth = write_to_file(CT, 'cell_count_6_receptors.csv')
print("Cell growth by day with 6 receptor clusters = ", cell_growth)


cTT, CT = cell_count(1/(24*3.5) , [24], e0_lst)
cell_growth = write_to_file(CT, 'cell_count_24_receptors.csv')
print("Cell growth by day with 24 receptor clusters = ", cell_growth)