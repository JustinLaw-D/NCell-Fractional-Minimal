import sys
sys.path.append('./../')

from NCell import NCell
import numpy as np
R = 6371 # radius of earth in km
dh = 50 # height of band (km)
alts = np.arange(600, 900+dh/2, dh)
alt_edges = np.arange(575, 925+1, dh)
V = 4*np.pi*dh*(R+alts)**2 # volume of band
S_i = []*len(alts)
D_i = []*len(alts)
lam = []
N_i = np.zeros(len(alts))
for i in range(len(alts)):
    S_i.append([0,0,0])
    D_i.append([0,0,0])
    N_i[i] = 2.5e-8*V[i]
    lam.append([0,0,0])
lam[0][0] = 1000
lam[3][1] = 500
lam[-1][2] = 100
T = 60
atmosphere = NCell(S_i, D_i, N_i, alt_edges, lam)

atmosphere.run_sim_precor(T, dt_min=1/10000, upper=True)
atmosphere.save("./", "test_save_NCellSystem", gap=0.1, force=True)
