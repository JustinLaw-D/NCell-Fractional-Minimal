# Test of stable single-cell system, using drag lifetime from JASON model

import sys
sys.path.append('./../')

from NCell import NCell
from Events import Event, ExplEvent
import numpy as np

class MyEvent(Event):

    def run_event(self, S, S_d, D, R, N, logL_edges, chi_edges):
        dN = np.zeros((len(logL_edges)-1, len(chi_edges)-1))
        dN[0,0] = 1e6
        return 0, 0, 0, 0, dN, [],[]

R = 6371 # radius of earth in km
alt = 600 # altitude of Starlink satellites (km)
dh = 25 # height of band (km)
V = 4*np.pi*dh*(R+alt)**2 # volume of band
S_i = [0]
S_di = [0]
D_i = [0]
N_i = 2.5e-8*V
lam = 2000
T = 50
events = [MyEvent(alt, time=[1,10])]#, ExplEvent(alt, [(1,'rb',10)], freq=1)]
atmosphere = NCell([S_i], [S_di], [D_i], [N_i], [alt], [587.5,612.5], [lam], tau_do=[[2]], events=events)
atmosphere.run_sim_precor(T, dt_min=1/10000)
atmosphere.save("./", "test_save_Discrete", gap=0.01)
