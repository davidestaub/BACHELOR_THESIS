#!/bin/python
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import rcParams
from matplotlib import ticker
rcParams.update({'figure.autolayout': True, 'font.size':8})
import itertools

import re
import pandas as pd
from os import listdir
import os
import glob
#os.chdir()

data = pd.read_csv("../../data.csv")

marker = itertools.cycle(('+','.','o','*'))
line = itertools.cycle((':','-.','--'))


fig_V, ax_V = plt.subplots()
plt.tight_layout()

ax_V.plot(data['time step'], data['Total Energy'], marker=next(marker),linestyle=next(line),label='Total Energy')
ax_V.set_xlabel("time step",fontdict=None,labelpad=0)
ax_V.set_ylabel("total energy",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Total Energy ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('total_energy.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
ax_V.plot(data['time step'], data['Total Energy'], marker=next(marker),linestyle=next(line),label='Total Energy')
ax_V.plot(data['time step'], data['Kinetic Energy'], marker=next(marker),linestyle=next(line),label='Kinetic Energy')
ax_V.plot(data['time step'], data['Potential Energy'], marker=next(marker),linestyle=next(line),label='Potential Energy')
ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
ax_V.set_ylabel("Energy [J]",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Energies ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('energy.pdf',format='pdf')


fig_V, ax_V = plt.subplots()
ax_V.plot(data['time step'], data['Velocity X'], marker=next(marker),linestyle=next(line),label='Velocity X')
ax_V.plot(data['time step'], data['Velocity Y'], marker=next(marker),linestyle=next(line),label='Velocity Y')
ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
ax_V.set_ylabel("Velocity ",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocities ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('velocities.pdf',format='pdf')

