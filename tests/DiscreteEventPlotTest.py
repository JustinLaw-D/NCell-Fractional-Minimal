import sys
sys.path.append("./../")

import numpy as np
from NCell import *

atmosphere = NCell.load("test_save_Discrete/")

t = atmosphere.get_t()
S = atmosphere.get_S()[0][0]
S_d = atmosphere.get_SD()[0][0]
D = atmosphere.get_D()[0][0]
N = atmosphere.get_N()[0]
C = atmosphere.get_C()[0]
T = t[-1]

import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()
ax1.set_xlabel('time (yr)')
ax1.set_ylabel('log(number)')
ax1.set_yscale('log')
ax1.plot(t, S, label='S')
ax1.plot(t, S_d, label='S_d')
ax1.plot(t, D, label='D')
ax1.plot(t, N, label='N')
ax1.set_ylim(100, 1e9)
ax1.set_xlim(0,T)
ax1.legend()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('log(collisions)')  # we already handled the x-label with ax1
ax2.plot(t, C, label='C', color='k')
ax2.set_ylim(1e-3, 1e8)
ax2.set_yscale('log')
ax2.legend()

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

