#!/bin/python
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import rcParams
from matplotlib import ticker
rcParams.update({'figure.autolayout': True, 'font.size':8})
import itertools
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import re
import pandas as pd
from os import listdir
import os
import glob
#os.chdir()

data = pd.read_csv("../../data.csv")
data2 = pd.read_csv("../../data2.csv")

velocities = pd.read_csv("../../velocities.csv")
print(velocities)

marker = itertools.cycle(('+','.','o','*'))
line = itertools.cycle((':','-.','--'))

fig_V, ax_V = plt.subplots()
plt.tight_layout()

ax_V.plot(velocities['time_step'], velocities['VELX'], marker=next(marker),linestyle=next(line),label='Velocity in x')
ax_V.plot(velocities['time_step'], velocities['VELY'], marker=next(marker),linestyle=next(line),label='Velocity in y')
ax_V.plot(velocities['time_step'], velocities['friction'], marker=next(marker),linestyle=next(line),label='friction in x')
ax_V.set_xlabel("time step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Velocity",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('Velocity.pdf',format='pdf')


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
ax_V.plot(data['time step'], data['E'], marker=next(marker),linestyle=next(line),label='E')
ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
ax_V.set_ylabel("E ",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("E ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('E.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
ax_V.plot(data['time step'], data['SPRING FORCE X'], marker=next(marker),linestyle=next(line),label='SPRING FORCE X')
#ax_V.plot(data['time step'], data['SPRING FORCE Y'], marker=next(marker),linestyle=next(line),label='SPRING FORCE Y')
ax_V.plot(data['time step'], data['DISTANCE'], marker=next(marker),linestyle=next(line),label='DISTANCE')
ax_V.plot(data['time step'], data['DERIVATIVE'], marker=next(marker),linestyle=next(line),label='DERIVATIVE')
ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
ax_V.set_ylabel("SPRING FORCE ",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("SPRING FORCE ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('SPRINGFORCE.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
ax_V.plot(data['time step'], data['X'], marker=next(marker),linestyle=next(line),label=' X')
ax_V.plot(data['time step'], data['Y'], marker=next(marker),linestyle=next(line),label='Y')
ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
ax_V.set_ylabel("Position",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Position ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('Position.pdf',format='pdf')

fig = plt.figure()
ax = Axes3D(fig)
print(data2)
# Plot the surface.
X = data2["x"]
Y = data2["y"]
Z1 = data2["z1"]
Z2 = data2["z2"]
Z = pd.DataFrame()
Z.append(Z1)
Z.append(Z2)



#X,Y = np.meshgrid(X,Y)
#surf = ax.plot_surface(X,Y,Z1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.set_xlabel("x ",fontdict=None,labelpad=0)
#ax.set_ylabel("y",fontdict=None,labelpad=-5, rotation=0,)
#ax.set_xlim(-400,400)
#ax.set_ylim(-400,400)
#ax.set_zlim(-10000, 10000)
#ax.set_title("Newton ", fontdict= None, loc = 'center',pad=5)
#fig.savefig('Newton.pdf',format='pdf')

#data3 = pd.read_csv("../../data3.csv")
#print(data3)



#x_best,y_best,z_best = np.loadtxt(open("../../best.csv"), delimiter=",", skiprows=1,unpack=True)
#print(x_best)
#print(y_best)
#print(z_best)
#best_data = pd.read_csv("../../best.csv")
#print(best_data)

#x,y,z = np.loadtxt(open("../../data3.csv"), delimiter=",", skiprows=1,unpack=True)
#fig = plt.figure()
#ax = Axes3D(fig)
#surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.1)
#ax.scatter(x_best,y_best,z_best)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.savefig('teste.pdf')



#dataFD = pd.read_csv("../../FD.csv")
#j = dataFD[(dataFD.i == 0 )].index
#print(dataFD)
#dataFD.drop(j)

#fig_V, ax_V = plt.subplots()
#ax_V.plot(dataFD['i'], dataFD['F_X'], marker=next(marker),linestyle=next(line),label='F X')
#ax_V.plot(dataFD['i'], dataFD['DEDX'], marker=next(marker),linestyle=next(line),label='E dx with FD')
#ax_V.plot(dataFD['i'], dataFD['Difference'], marker=next(marker),linestyle=next(line),label='Difference')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("FD ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('FD.pdf',format='pdf')




#fig_V, ax_V = plt.subplots()
#ax_V.plot(dataFD['i'], dataFD['FORCE_JACOBIAN'], marker=next(marker),linestyle=next(line),label='Force JAcobian')
#ax_V.plot(dataFD['i'], dataFD['DFDX'], marker=next(marker),linestyle=next(line),label='F dx with FD')
#ax_V.plot(dataFD['i'], dataFD['DIFFERENCE2'], marker=next(marker),linestyle=next(line),label='Difference')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("FD2 ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('FD2.pdf',format='pdf')

#fig_V, ax_V = plt.subplots()
#ax_V.plot(dataFD['i'], dataFD['FORCE_JACOBIAN'], marker=next(marker),linestyle=next(line),label='Force JAcobian')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("FD3 ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('FD3.pdf',format='pdf')

#fig_V, ax_V = plt.subplots()
#ax_V.plot(dataFD['i'], dataFD['DFDX'], marker=next(marker),linestyle=next(line),label='F dx with FD')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("FD4 ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('FD4.pdf',format='pdf')

#test_e = pd.read_csv("../../test_e.csv")

#print(test_e)
#x,y,z = np.loadtxt(open("../../test_e.csv"), delimiter=",", skiprows=1,unpack=True)
#fig = plt.figure()
#ax = Axes3D(fig)
#surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.1)
#ax.scatter(x_best,y_best,z_best)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.savefig('test_e.pdf')


plt.show()
