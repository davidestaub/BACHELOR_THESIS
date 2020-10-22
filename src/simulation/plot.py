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

position = pd.read_csv("../../position.csv")
rotator_position = pd.read_csv("../../rotator_position.csv")
triangle_position = pd.read_csv("../../triangle_position.csv")
print("triangle positions")
print(triangle_position)
print("rotator_position")
print(rotator_position)
FMA = pd.read_csv("../../FMA.csv")
brownian = pd.read_csv("../../brownian.csv")
drag = pd.read_csv("../../drag.csv")
print("brownian")
print(brownian)
#print(position)
ehd = pd.read_csv("../../ehd.csv")
spring = pd.read_csv("../../spring.csv")
print("spring")
print(spring)
print("ehd")
print(ehd)
amount = pd.read_csv("../../amount.csv")
print(amount)
print("position")
print(position)

data = pd.read_csv("../../data.csv")
data2 = pd.read_csv("../../data2.csv")
v_w_data = pd.read_csv("../../v_w.csv")

print("v_W_data")
pd.set_option("display.max_rows", None, "display.max_columns", None)
print(v_w_data)

v_E_data = pd.read_csv("../../v_E.csv")
#print(v_E_data)

velocities = pd.read_csv("../../velocities.csv")
#print(velocities)

marker = itertools.cycle(('+','.','o','*'))
line = itertools.cycle((':','-.','--'))


fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(position['Position in X']* pow(10,6), position['Position in Y']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.plot(position['Position in X 2'] * pow(10,6), position['Position in Y 2']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
ax_V.plot(position['Position in X 3']* pow(10,6), position['Position in Y 3']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
#ax_V.plot(position['Position in X 4']* pow(10,6), position['Position in Y 4']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))

ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Dumbbell Trajectories ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('path.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(rotator_position['Position in X']* pow(10,6), rotator_position['Position in Y']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.plot(rotator_position['Position in X 2'] * pow(10,6), rotator_position['Position in Y 2']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)

ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Rotator Trajectories ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('path.pdf',format='pdf')

cm = plt.cm.get_cmap('RdYlBu')
fig_V, ax_V = plt.subplots()
plt.tight_layout

sc = plt.scatter(x=position['Position in X'],y=position['Position in Y'],c=position['Time Step'], cmap=cm)
plt.colorbar(sc)
plt.xlim(-0.0081,0.001)
plt.ylim(-0.0008,0.00025)
plt.xlabel("X",fontdict=None,labelpad=0)
plt.ylabel("Y",fontdict=None,labelpad=-5, rotation=0,)
plt.legend(loc='best',fontsize=8)
plt.title("path ", fontdict= None, loc = 'center',pad=5)
plt.savefig('path_with_colors.pdf',format='pdf')
#ax_V.scatter(position['Position in X'], position['Position in Y'], c = position['Time Step'])
#ax_V.set_xlim(-0.0012,0.0006)
#ax_V.set_ylim(0.00,0.0012)
#fig.colorbar()
#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))

#ax_V.set_xlabel("X",fontdict=None,labelpad=0)
#ax_V.set_ylabel("Y",fontdict=None,labelpad=-5, rotation=0,)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("path ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('path_with_colors.pdf',format='pdf')


fig_V, ax_V = plt.subplots()
plt.tight_layout
sc = plt.scatter(x=rotator_position['Position in X']* pow(10,6),y=rotator_position['Position in Y']* pow(10,6),c=rotator_position['Time Step'] * 2.0 * pow(10,-3), cmap=cm,s = 0.8)
colorb = plt.colorbar(sc)
colorb.ax.get_yaxis().labelpad = 15
colorb.ax.set_ylabel('Time [s]', rotation=270)
plt.xlabel("x(t) [µm]",fontdict=None,labelpad=0)
plt.ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
plt.ticklabel_format(axis='both', style = 'sci')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Rotator Trajectories", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('rot_path.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(position['Time Step'], position['Position in X'], marker=next(marker),linestyle=next(line))
#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))

ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Position",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Position new ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('pos_new.pdf',format='pdf')


fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(triangle_position['Triangle Position in X']* pow(10,6), triangle_position['Triangle Position in Y']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.plot(triangle_position['Triangle Position in X 2'] * pow(10,6), triangle_position['Triangle Position in Y 2']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
ax_V.plot(triangle_position['Triangle Position in X 3'] * pow(10,6), triangle_position['Triangle Position in Y 3']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
ax_V.plot(triangle_position['COM Position in X 4'] * pow(10,6), triangle_position['COM Position in Y 4']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Triangle Trajectories ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('tri_path.pdf',format='pdf')


fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(FMA['time step'], FMA['F-M*a'], marker=next(marker),linestyle=next(line))
#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))

ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("F-M*a",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("F = M * A ? ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('FMA.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(brownian['time step'], brownian['Brownian Force X'], marker=next(marker),linestyle=next(line))
ax_V.plot(brownian['time step'], brownian['Brownian Force Y'], marker=next(marker),linestyle=next(line))

#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Brownian Force",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Brownian Force over time ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('brownian.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(amount['time step'], amount['Amount X'], marker=next(marker),linestyle=next(line))
ax_V.plot(amount['time step'], amount['Amount Y'], marker=next(marker),linestyle=next(line))
ax_V.plot(amount['time step'], amount['Amount T'], marker=next(marker),linestyle=next(line))

#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Amount",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Direction", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('dir.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(amount['time step'], amount['Norm'], marker=next(marker),linestyle=next(line))
#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Direction Norm",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Direction Norm", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('dir_norm.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(drag['time step'], drag['Drag X'], marker=next(marker),linestyle=next(line))
ax_V.plot(drag['time step'], drag['Drag Y'], marker=next(marker),linestyle=next(line))
ax_V.plot(drag['time step'], drag['Drag Px'], marker=next(marker),linestyle=next(line))
ax_V.plot(drag['time step'], drag['Drag Py'], marker=next(marker),linestyle=next(line))
#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Drag",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Drag", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('drag.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(amount['Amount X'], amount['Amount Y'], marker=next(marker),linestyle=next(line))


#ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Amount X",fontdict=None,labelpad=0)
ax_V.set_ylabel("Amount Y",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Direction", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('dir2.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(ehd['time step'], ehd['EHD Force X'], marker=next(marker),linestyle=next(line))
ax_V.plot(ehd['time step'], ehd['EHD Force Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("EHD Force",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("EHD Force over time ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('ehd.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(spring['time step'], spring['Spring Force X'], marker=next(marker),linestyle=next(line))
ax_V.plot(spring['time step'], spring['Spring Force Y'], marker=next(marker),linestyle=next(line))
ax_V.set_xlabel("Time Step",fontdict=None,labelpad=0)
ax_V.set_ylabel("Spring Force",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Spring Force over time ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('spring.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot((v_w_data['Frequency']), v_w_data[' Velocity'] *-1 * pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.set_xscale('log')
ax_V.axhline(linewidth=4, color='k')
ax_V.set_xlabel("Frequency",fontdict=None,labelpad=0)
ax_V.set_ylabel("Velocity",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity vs Frequency in micrometer/s ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('Velocity_frequency.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(v_E_data['E Field'], v_E_data[' Velocity'] , marker=next(marker),linestyle=next(line))

ax_V.set_xlabel("E Field",fontdict=None,labelpad=0)
ax_V.set_ylabel("Velocity",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity vs E Fieldin micrometer/s ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('Velocity_E_Field.pdf',format='pdf')



#fig_V, ax_V = plt.subplots()
#plt.tight_layout()

#ax_V.plot(velocities['time_step'], velocities['VELX'], marker=next(marker),linestyle=next(line),label='Velocity in x')
#ax_V.plot(velocities['time_step'], velocities['VELY'], marker=next(marker),linestyle=next(line),label='Velocity in y')
#ax_V.plot(velocities['time_step'], velocities['friction'], marker=next(marker),linestyle=next(line),label='friction in x')
#ax_V.set_xlabel("time step",fontdict=None,labelpad=0)
#ax_V.set_ylabel("Velocity",fontdict=None,labelpad=-5, rotation=0,)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("Velocity ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('Velocity.pdf',format='pdf')


#fig_V, ax_V = plt.subplots()
#plt.tight_layout()

#ax_V.plot(data['time step'], data['Total Energy'], marker=next(marker),linestyle=next(line),label='Total Energy')
#ax_V.set_xlabel("time step",fontdict=None,labelpad=0)
#ax_V.set_ylabel("total energy",fontdict=None,labelpad=-5, rotation=0,)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("Total Energy ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('total_energy.pdf',format='pdf')

#fig_V, ax_V = plt.subplots()
#ax_V.plot(data['time step'], data['Total Energy'], marker=next(marker),linestyle=next(line),label='Total Energy')
#ax_V.plot(data['time step'], data['Kinetic Energy'], marker=next(marker),linestyle=next(line),label='Kinetic Energy')
#ax_V.plot(data['time step'], data['Potential Energy'], marker=next(marker),linestyle=next(line),label='Potential Energy')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.set_ylabel("Energy [J]",fontdict=None,labelpad=-5, rotation=0,)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("Energies ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('energy.pdf',format='pdf')




#fig_V, ax_V = plt.subplots()
#ax_V.plot(data['time step'], data['E'], marker=next(marker),linestyle=next(line),label='E')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.set_ylabel("E ",fontdict=None,labelpad=-5, rotation=0,)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("E ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('E.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
ax_V.plot(data['time step'], data['SPRING FORCE X'], marker=next(marker),linestyle=next(line),label='SPRING FORCE X')
ax_V.plot(data['time step'], data['SPRING FORCE Y'], marker=next(marker),linestyle=next(line),label='SPRING FORCE Y')
#ax_V.plot(data['time step'], data['DISTANCE'], marker=next(marker),linestyle=next(line),label='DISTANCE')
#ax_V.plot(data['time step'], data['DERIVATIVE'], marker=next(marker),linestyle=next(line),label='DERIVATIVE')
ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
ax_V.set_ylabel("SPRING FORCE ",fontdict=None,labelpad=-5, rotation=0,)
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("SPRING FORCE ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('SPRINGFORCE.pdf',format='pdf')

#fig_V, ax_V = plt.subplots()
#ax_V.plot(data['time step'], data['X'], marker=next(marker),linestyle=next(line),label=' X')
#ax_V.plot(data['time step'], data['Y'], marker=next(marker),linestyle=next(line),label='Y')
#ax_V.set_xlabel("time step ",fontdict=None,labelpad=0)
#ax_V.set_ylabel("Position",fontdict=None,labelpad=-5, rotation=0,)
#ax_V.legend(loc='best',fontsize=8)
#ax_V.grid(True,linestyle='-')
#ax_V.set_title("Position ", fontdict= None, loc = 'center',pad=5)
#fig_V.savefig('Position.pdf',format='pdf')

#fig = plt.figure()
#ax = Axes3D(fig)
#print(data2)
# Plot the surface.
#X = data2["x"]
#Y = data2["y"]
#Z1 = data2["z1"]
#Z2 = data2["z2"]
#Z = pd.DataFrame()
#Z.append(Z1)
#Z.append(Z2)














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
