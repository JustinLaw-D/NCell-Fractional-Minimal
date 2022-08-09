import sys
sys.path.append('./../')

from NCell import NCell
import numpy as np

R = 6371 # radius of earth in km
dh = 50 # height of band (km)
alts = np.arange(600, 910, dh)
alt_edges = np.arange(575, 925+1, dh)
V = 4*np.pi*dh*(R+alts)**2 # volume of band
S_i=[]
S_di=[]
D_i=[]
tau_do=[]
for i in range(len(alts)):
    S_i.append([0])
    S_di.append([0])
    D_i.append([0])
    tau_do.append([10])
N_i = np.zeros(len(alts), dtype=np.int64)
S_di[-1][0] = 100
target_alts = [500]
lam = [0]
T = 30
atmosphere = NCell(S_i, S_di, D_i, N_i, target_alts, alt_edges, lam, tau_do=tau_do)
atmosphere.run_sim_precor(T, dt_min=1/1000, upper=False)
atmosphere.save('./', "test_save_NDeorbitFlow", gap=0.05)
