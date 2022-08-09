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
S_di=[]
D_i=[]
up_time=[]
for i in range(len(alts)):
    S_i.append([0])
    S_di.append([0])
    D_i.append([0])
    up_time.append([0.5])
N_i = np.zeros(len(alts), dtype=np.int64)
target_alts = [alts[-1]]
lam = [50]
T = 100
atmosphere = NCell(S_i, S_di, D_i, N_i, target_alts, alt_edges, lam, up_time=up_time)

atmosphere.run_sim_precor(T, upper=False)
atmosphere.save('./', "test_save_NLaunchFlow", gap=0.05, force=True)
