import sys
sys.path.append('./../')

from NCell import NCell
import numpy as np
R = 6371 # radius of earth in km
alt1, alt2 = 600, 625 # altitude of Starlink satellites (km)
dh = 25 # height of bands (km)
alt_edges = [587.5, 612.5, 637.5]
V1, V2 = 4*np.pi*dh*(R+alt1)**2, 4*np.pi*dh*(R+alt2)**2 # volume of bands
S_i = [[0,0,0],[0,0,0]]
S_di = [[0,0,0],[0,0,0]]
D_i = [[0,0,0],[0,0,0]]
N_i1, N_i2 = int(0), int(0)
lam = [0,0,0]
target_alts = [alt2,alt2,alt2]
R_i = [[0,0,0],[50,1e2,2e2]]
T = 50
atmosphere = NCell(S_i, S_di, D_i, [N_i1, N_i2], target_alts, alt_edges, lam, R_i=R_i)

atmosphere.run_sim_precor(T, upper=False)
atmosphere.save('./', 'test_save_2CellNTypeRocket', gap=0.01)
