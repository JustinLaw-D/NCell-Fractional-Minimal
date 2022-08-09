import sys
sys.path.append('./../')

from NCell import NCell
import numpy as np

atmosphere = NCell.load('./test_save_NCellSystem/')
t = atmosphere.get_t()
T = t[-1]
S_data = atmosphere.get_S()
Sd_data = atmosphere.get_SD()
D_data = atmosphere.get_D()
N = atmosphere.get_N()
C = atmosphere.get_C()
S = []
S_d = []
D = []
for i in range(atmosphere.num_cells):
    S_loc, Sd_loc, D_loc = [], [], []
    for j in range(len(t)):
        S_loc.append(0)
        Sd_loc.append(0)
        D_loc.append(0)
        for k in range(atmosphere.num_sat_types):
            S_loc[j] += S_data[i][k][j]
            Sd_loc[j] += Sd_data[i][k][j]
            D_loc[j] += D_data[i][k][j]
    S.append(S_loc)
    S_d.append(Sd_loc)
    D.append(D_loc)

import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()
ax1.set_xlabel('time (yr)')
ax1.set_ylabel('log(number)')
ax1.set_yscale('log')
for i in range(atmosphere.num_cells):
    #ax1.plot(t, S[i], label='S'+str(i))
    #ax1.plot(t, S_d[i], label='SD'+str(i))
    #ax1.plot(t, D[i], label='D'+str(i))
    #ax1.plot(t, N[i], label='N'+str(i))
    ax1.plot(t, C[i], label='C'+str(i))
ax1.set_ylim(1, 1e9)
ax1.set_xlim(0,T)
ax1.legend()

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
