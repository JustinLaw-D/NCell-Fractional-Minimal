import sys
sys.path.append('./../')
sys.path.append('./../catalog_data')

from NCell import NCell
from data_utilities import *
import numpy as np

R = 6371 # radius of earth in km
dh = 50 # height of band (km)
alt_bins = np.arange(600-dh/2, 900+dh/2+1, dh)
N_i = []
S_i=[]
D_i=[]
lam=[]
for i in range(len(alt_bins)-1):
    S_i.append([0])
    D_i.append([0])
    lam.append([0])
    N_i.append(0)
N_i[-1] = 1000
T = 100
atmosphere = NCell(S_i, D_i, N_i, alt_bins, lam)

atmosphere.run_sim_precor(T, upper=False)
atmosphere.save('./', "test_save_NDebrisFlow", gap=0.05, force=True)
