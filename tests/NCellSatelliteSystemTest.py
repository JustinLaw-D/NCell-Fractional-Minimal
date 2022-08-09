import sys
sys.path.append('./../')

from NCell import NCell
import numpy as np
R = 6371 # radius of earth in km
dh = 50 # height of band (km)
alts = np.arange(600, 900+dh/2, dh)
alt_edges = np.arange(575, 925+1, dh)
V = 4*np.pi*dh*(R+alts)**2 # volume of band
S_i = [[0,0,0]]*len(alts)
D_i = [[0,0,0]]*len(alts)
S_di = [[0,0,0]]*len(alts)
N_i = np.zeros(len(alts))
for i in range(len(alts)):
    N_i[i] = 2.5e-8*V[i]
lam = [1000, 1000, 1000]
T = 60
atmosphere = NCell(S_i, S_di, D_i, N_i, [600, 750, 900],  alt_edges, lam)

atmosphere.run_sim_precor(T, dt_min=1/1000, upper=True)
atmosphere.save("./", "test_save_NCellSystem", gap=0.1, force=True)
