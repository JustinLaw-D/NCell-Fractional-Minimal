import sys
sys.path.append('./../')

from NCell import NCell
import numpy as np
R = 6371 # radius of earth in km
dh = 50 # height of band (km)
alt_edges = np.arange(575, 926, dh)
alts = np.arange(600, 910, dh)
V = 4*np.pi*dh*(R+alts)**2 # volume of band
S_i=[]
D_i=[]
lam=[]
for i in range(len(alts)):
    S_i.append([0])
    D_i.append([0])
    lam.append([0])
N_i = np.zeros(len(alts), dtype=np.int64)
lam[-1][0], lam[0][0] = 200, 100
T = 100
atmosphere = NCell(S_i, D_i, N_i, alt_edges, lam)

atmosphere.run_sim_euler(T, dt=0.001, upper=False)
atmosphere.save('./', "test_save_NLaunchFlow", gap=0.05, force=True)
