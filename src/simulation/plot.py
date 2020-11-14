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
from math import sqrt

import re
import pandas as pd
from os import listdir
import os
import glob
#os.chdir()

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')


# Function to find the circle on
# which the given three points lie
def findCircle(x1, y1, x2, y2, x3, y3) :
    x12 = x1 - x2;
    x13 = x1 - x3;

    y12 = y1 - y2;
    y13 = y1 - y3;

    y31 = y3 - y1;
    y21 = y2 - y1;

    x31 = x3 - x1;
    x21 = x2 - x1;

    # x1^2 - x3^2
    sx13 = pow(x1, 2) - pow(x3, 2);

    # y1^2 - y3^2
    sy13 = pow(y1, 2) - pow(y3, 2);

    sx21 = pow(x2, 2) - pow(x1, 2);
    sy21 = pow(y2, 2) - pow(y1, 2);

    f = (((sx13) * (x12) + (sy13) *
          (x12) + (sx21) * (x13) +
          (sy21) * (x13)) // (2 *
                              ((y31) * (x12) - (y21) * (x13))));

    g = (((sx13) * (y12) + (sy13) * (y12) +
          (sx21) * (y13) + (sy21) * (y13)) //
         (2 * ((x31) * (y12) - (x21) * (y13))));

    c = (-pow(x1, 2) - pow(y1, 2) -
         2 * g * x1 - 2 * f * y1);

    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
    # where centre is (h = -g, k = -f) and
    # radius r as r^2 = h^2 + k^2 - c
    h = -g;
    k = -f;
    sqr_of_r = h * h + k * k - c;

    # r is the radius
    r = round(sqrt(sqr_of_r), 10);
   # print("Centre = (", h, ", ", k, ")");
   # print("Radius = ", r);
    return r






marker = itertools.cycle(('+','.','o','*'))
line = itertools.cycle((':','-.','--'))

sigma_m_data= pd.read_csv('../../sigma_m.csv')
fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(sigma_m_data['Medium conductivity'], sigma_m_data['Velocity']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.set_xlabel("medium conductivity",fontdict=None,labelpad=0)
ax_V.set_ylabel("velocity",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity vs. Medium conductivity", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('sigma_m.pdf',format='pdf')



sigma_p_data= pd.read_csv('../../sigma_p.csv')
fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(sigma_p_data['Particle conductivity'], sigma_p_data['Velocity']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.set_xlabel("particle conductivity",fontdict=None,labelpad=0)
ax_V.set_ylabel("velocity",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity vs. Particle conductivity", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('sigma_p.pdf',format='pdf')


epsilon_p_data= pd.read_csv('../../epsilon_p.csv')
fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(epsilon_p_data['Relative particle permittivity'], epsilon_p_data['Velocity']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.set_xlabel("particle permittivity",fontdict=None,labelpad=0)
ax_V.set_ylabel("velocity",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity vs. Particle permittivity", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('epsilon_p.pdf',format='pdf')


epsilon_m_data= pd.read_csv('../../epsilon_m.csv')
fig_V, ax_V = plt.subplots()
plt.tight_layout
ax_V.plot(epsilon_m_data['Relative medium permittivity'], epsilon_m_data['Velocity']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.set_xlabel("Medium permittivity",fontdict=None,labelpad=0)
ax_V.set_ylabel("velocity",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
ax_V.grid(True,linestyle='-')
ax_V.set_title("Velocity vs. medium permittivity", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('epsilon_m.pdf',format='pdf')



##begin: plot curvature of experimental dumbell data####

column_names_velocity = ["Velocity"]
velocities = pd.DataFrame(columns= column_names_velocity)

experimental_data_dumbbell = pd.read_csv('../../Active_dumbbels_trajectories.txt',sep="   ", header=None,skiprows=1)
experimental_data_dumbbell.columns = ['x','y','frame_id','particle_id']
experimental_data_dumbbell['x'] = experimental_data_dumbbell['x'] * 0.41
experimental_data_dumbbell['y'] = experimental_data_dumbbell['y'] * 0.41
#print("hello")
#print(experimental_data_dumbbell)
experimental_data_dumbbell['Time passed [s]'] = (experimental_data_dumbbell['frame_id'] - 1) * 100 * pow(10,-3)
fig_V, ax_V = plt.subplots()
plt.tight_layout
column_names = ["radius","curvature","velocity"]
curvature_data = pd.DataFrame(columns = column_names)
unique_ids = experimental_data_dumbbell.particle_id.unique()
print(unique_ids)
n = unique_ids.shape[0]
n=15
print(n)
time = [100,1000]
print(time)
for k in time:
    for i in range(1,15,1):
        print(i)
        print(i)
        is_current_id = experimental_data_dumbbell['particle_id'] == i
        current_experimental_data = experimental_data_dumbbell[is_current_id]
        current_experimental_data.reset_index(inplace=True)
        column_names_velocity = ["Velocity"]
        velocities = pd.DataFrame(columns= column_names_velocity)


        fig_C, ax_C = plt.subplots()
        ax_C.plot(current_experimental_data['x'], current_experimental_data['y'], marker=next(marker),linestyle=next(line),markersize = 0.8)
        ax_C.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
        ax_C.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
        ax_C.ticklabel_format(axis='both', style = 'sci')
        ax_C.legend(loc='best',fontsize=8)
        ax_C.grid(True,linestyle='-')
        ax_C.set_title("Experimental Dumbbell Trajectories %i" %i + "time=%i" %k, fontdict= None, loc = 'center',pad=5)
        fig_C.savefig('plot_results/exp_path_%i' %i +'_time=%i' %k +'.pdf',format='pdf')

        cm = plt.cm.get_cmap('RdYlBu')
        fig_V, ax_V = plt.subplots()
        plt.tight_layout
        sc = plt.scatter(x=current_experimental_data['x'],y=current_experimental_data['y'],c=current_experimental_data['Time passed [s]'], cmap=cm, s=0.7)
        plt.colorbar(sc)
        plt.xlabel("x(t) [µm]",fontdict=None,labelpad=0,fontsize=10)
        plt.ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,fontsize=10)
        plt.legend(loc='best',fontsize=8)
        plt.title("Experimental Dumbbell Trajectories  %i" %i + "time=%i" %k, fontdict= None, loc = 'center',pad=5)
        plt.savefig('plot_results/Experimental Dumbbell Trajectories %i' %i +'_time=%i' %k +'.pdf',format='pdf')




    #

        for j in range(0,current_experimental_data.shape[0]-round(2 * k/100),): #loop over trajectory for one particle
            x1 =current_experimental_data.loc[j,'x']
            x2 =current_experimental_data.loc[j+round(1 * k/100),'x']
            x3 =current_experimental_data.loc[j+round(2 * k/100),'x']
            y1 =current_experimental_data.loc[j,'y']
            y2 =current_experimental_data.loc[j+round(1 * k/100),'y']
            y3 =current_experimental_data.loc[j+round(2 * k/100),'y']

       #     print(" current exp")
      #      print(current_experimental_data.loc[j])
     #       print(" next exp")
            #print(current_experimental_data.loc[j+1])
            #print(" nextnext exp")
           # print(current_experimental_data.loc[j+2])
            r = findCircle(x1,y1,x2,y2,x3,y3)
            c = 1.0/r
            #print(" r")

            t1 =pow(x2-x1,2)+pow(y2-y1,2)
            t2 =pow(x3-x2,2)+pow(y3-y2,2)
            vel1 = sqrt(t1)/(k * pow(10,-3))
            vel2 = sqrt(t2)/(k * pow(10,-3))
            vel = (vel1 + vel2) * 0.5
            pdvelcurrent1 = pd.DataFrame([[vel1]],columns = column_names_velocity)
            pdvelcurrent2 = pd.DataFrame([[vel2]],columns = column_names_velocity)
            ptmp = velocities.append(pdvelcurrent1,ignore_index=True)
            velocities = ptmp
            ptmp = velocities.append(pdvelcurrent2,ignore_index=True)
            velocities = ptmp
            pd2 = pd.DataFrame([[r,c,vel]], columns = column_names)
            pd3 = curvature_data.append(pd2,ignore_index=True)
            curvature_data = pd3
#
#
#
#
#
        fig_V, ax_V = plt.subplots()
        plt.tight_layout
        num_bins = 100
        n, bins, patches = ax_V.hist(velocities["Velocity"], num_bins, density=1)
        sigma = velocities["Velocity"].std()
        mu = velocities["Velocity"].mean()
        print("mu")
        print(mu)
        print("sigma")
        print(sigma)
        y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
             np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
        ax_V.plot(bins, y, '--')
        ax_V.set_xlabel("velocity [µm/s]",fontdict=None,labelpad=0, rotation=0,fontsize=10)
        ax_V.set_ylabel("Probability Density",fontdict=None,labelpad=0,fontsize=10)
        ax_V.ticklabel_format(axis='both', style = 'sci')
        #ax_V.legend(loc='best',fontsize=8)
        ax_V.grid(True,linestyle='-')
        ax_V.set_title('Histogram for Dumbbell %i' %i + ' and time %i' %k + 'mean is %f' %mu + 'and std is %f' %sigma,fontdict= None, loc = 'center',pad=5)
        #ax_V.set_title("Experimental Curvature vs. Velocity", fontdict= None, loc = 'center',pad=5)
        fig_V.tight_layout()
        fig_V.savefig('plot_results/exp_velocity_histogram %i' %i +' and time %i.pdf' %k,format='pdf')
#
#
#
#         ax_V.plot(current_experimental_data['x'], current_experimental_data['y'], marker=next(marker),linestyle=next(line),markersize = 0.8)
#         ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
#         ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
#         ax_V.ticklabel_format(axis='both', style = 'sci')
#         ax_V.legend(loc='best',fontsize=8)
#         ax_V.grid(True,linestyle='-')
#         ax_V.set_title("Experimental Dumbbell Trajectories ", fontdict= None, loc = 'center',pad=5)
#         fig_V.savefig('exp_path.pdf',format='pdf')

#plt.hist(velocities["velocities"],100,density=True)

    fig_V, ax_V = plt.subplots()
    plt.tight_layout
    num_bins = 100
    n, bins, patches = ax_V.hist(velocities["Velocity"], num_bins, density=1)
    sigma = velocities["Velocity"].std()
    mu = velocities["Velocity"].mean()
    print("mu")
    print(mu)
    print("sigma")
    print(sigma)
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    ax_V.plot(bins, y, '--')
    ax_V.set_xlabel("velocity [µm/s]",fontdict=None,labelpad=0, rotation=0,fontsize=10)
    ax_V.set_ylabel("Probability Density",fontdict=None,labelpad=0,fontsize=10)
    ax_V.ticklabel_format(axis='both', style = 'sci')
    #ax_V.legend(loc='best',fontsize=8)
    ax_V.grid(True,linestyle='-')
    #ax_V.set_title("Experimental Curvature vs. Velocity", fontdict= None, loc = 'center',pad=5)
    fig_V.tight_layout()
    plt.show()
    fig_V.savefig('plot_results/exp_velocity_histogram_time_%i.pdf' %k,format='pdf')

### End: plot curvature of experimental dumbbell data ##


## Begin: plot curvature of simulated dumbbell data ##
column_names_velocity = ["Velocity"]
total_velocities = pd.DataFrame(columns= column_names_velocity)
simulated_data_dumbbell = pd.read_csv("../../curvature.csv")
simulated_data_dumbbell['x'] = simulated_data_dumbbell['x'] * pow(10,6)
simulated_data_dumbbell['y'] = simulated_data_dumbbell['y'] * pow(10,6)
fig_V, ax_V = plt.subplots()
plt.tight_layout
column_names = ["radius","curvature","velocity"]
curvature_data_sim = pd.DataFrame(columns = column_names)
unique_ids = simulated_data_dumbbell.dumbbell_index.unique()
n = unique_ids.shape[0]
time = [100,1000]
for k in time:
    for i in range(0,n,1):
        is_current_id = simulated_data_dumbbell['dumbbell_index'] == i
        current_simulated_data = simulated_data_dumbbell[is_current_id]
        current_simulated_data.reset_index(inplace=True)

        column_names_velocity = ["Velocity"]
        velocities = pd.DataFrame(columns= column_names_velocity)


        fig_C, ax_C = plt.subplots()
        ax_C.plot(current_simulated_data['x'], current_simulated_data['y'], marker=next(marker),linestyle=next(line),markersize = 0.8)
        ax_C.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
        ax_C.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
        ax_C.ticklabel_format(axis='both', style = 'sci')
        ax_C.legend(loc='best',fontsize=8)
        ax_C.grid(True,linestyle='-')
        ax_C.set_title("Simulated Dumbbell Trajectories ks=25 bm=0.3 beta=0.15 %i" %i + "time=%i" %k, fontdict= None, loc = 'center',pad=5)
        fig_C.savefig('plot_results/Simulated Dumbbell Trajectories ks=25 bm=0.3 beta=0.15_path_%i' %i +'_time=%i' %k +'.pdf',format='pdf')


        cm = plt.cm.get_cmap('RdYlBu')
        fig_V, ax_V = plt.subplots()
        plt.tight_layout
        sc = plt.scatter(x=current_simulated_data['x'],y=current_simulated_data['y'],c=current_simulated_data['time_passed'], cmap=cm, s=0.7)
        plt.colorbar(sc)
        plt.xlabel("x(t) [µm]",fontdict=None,labelpad=0,fontsize=10)
        plt.ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,fontsize=10)
        plt.legend(loc='best',fontsize=8)
        plt.title("Simulated Dumbbell Trajectories ks=25 bm=0.3 beta=0.15 %i" %i + "time=%i" %k, fontdict= None, loc = 'center',pad=5)
        plt.savefig('plot_results/Simulated Dumbbell Trajectories ks=25 bm=0.3 beta=0.15_path_%i' %i +'_time=%i' %k +'.pdf',format='pdf')


        for j in range(round(current_simulated_data.shape[0]/6),current_simulated_data.shape[0]-2,10): #loop over trajectory for one particle
            current_time_passed = current_simulated_data.loc[j,'time_passed']
            final_row= current_simulated_data.tail(1)
            max_time_passed = current_simulated_data['time_passed'].iloc[-1]
            if(max_time_passed - current_time_passed > (300 * pow(10,-3))):
                next_time_passed = current_time_passed + (100 * pow(10,-3))
                next_next_time_passed = next_time_passed + (100 * pow(10,-3))
                #valid_data_next = current_simulated_data[current_simulated_data['time_passed']==next_time_passed]
                valid_data_next = current_simulated_data.iloc[(current_simulated_data['time_passed']-next_time_passed).abs().argsort()[:1]]
                valid_data_next_next = current_simulated_data.iloc[(current_simulated_data['time_passed']-next_next_time_passed).abs().argsort()[:1]]
                #valid_data_next_next = current_simulated_data[current_simulated_data['time_passed']==next_next_time_passed]
                print("current")
                print(current_simulated_data.loc[j])
                print("next")
                print(valid_data_next)
                print("next next")
                print(valid_data_next_next)
                x1 =current_simulated_data.loc[j,'x']
                x2 =valid_data_next['x'].values
                x3 =valid_data_next_next['x'].values
                y1 =current_simulated_data.loc[j,'y']
                y2 =valid_data_next['y'].values
                y3 =valid_data_next_next['y'].values
                r = findCircle(x1,y1,x2,y2,x3,y3)
                print("r")
                print(r)
                print("                               ")
                c = 1.0/r
                t1 =pow(x2-x1,2)+pow(y2-y1,2)
                t2 =pow(x3-x2,2)+pow(y3-y2,2)
                vel1 = sqrt(t1)/(100 * pow(10,-3))
                vel2 = sqrt(t2)/(100 * pow(10,-3))
                vel = (vel1 + vel2) * 0.5
                pdvelcurrent1 = pd.DataFrame([[vel1]],columns = column_names_velocity)
                pdvelcurrent2 = pd.DataFrame([[vel2]],columns = column_names_velocity)
                ptmp = velocities.append(pdvelcurrent1,ignore_index=True)
                velocities = ptmp
                ptmp = velocities.append(pdvelcurrent2,ignore_index=True)
                velocities = ptmp
                pd2 = pd.DataFrame([[r,c,vel]], columns = column_names)
                pd3 = curvature_data_sim.append(pd2,ignore_index=True)
                curvature_data_sim = pd3
                total_velocities = total_velocities.append(velocities)
        fig_V, ax_V = plt.subplots()
        plt.tight_layout
        num_bins = 100
        n, bins, patches = ax_V.hist(velocities["Velocity"], num_bins, density=1)
        sigma = velocities["Velocity"].std()
        mu = velocities["Velocity"].mean()
        print("mu")
        print(mu)
        print("sigma")
        print(sigma)
        y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
        ax_V.plot(bins, y, '--')
        ax_V.set_xlabel("velocity [µm/s]",fontdict=None,labelpad=0, rotation=0,fontsize=10)
        ax_V.set_ylabel("Probability Density",fontdict=None,labelpad=0,fontsize=10)
        ax_V.ticklabel_format(axis='both', style = 'sci')
        #ax_V.legend(loc='best',fontsize=8)
        ax_V.grid(True,linestyle='-')
        ax_V.set_title('Histogram for Simulated Dumbbell ks=25 bm=0.3 beta=0.15 %i' %i + ' and time %i' %k + 'mean is %f' %mu + 'and std is %f' %sigma,fontdict= None, loc = 'center',pad=5)
        #ax_V.set_title("Experimental Curvature vs. Velocity", fontdict= None, loc = 'center',pad=5)
        fig_V.tight_layout()
        fig_V.savefig('plot_results/Simulated Dumbbell ks=25 bm=0.3 beta=0.15_velocity_histogram %i' %i +' and time %i.pdf' %k,format='pdf')





        ax_V.plot(current_simulated_data['x'], current_simulated_data['y'], marker=next(marker),linestyle=next(line),markersize = 0.8)
        ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
        ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
        ax_V.ticklabel_format(axis='both', style = 'sci')
        ax_V.grid(True,linestyle='-')
        ax_V.set_title("Simulated Dumbbell Trajectories ", fontdict= None, loc = 'center',pad=5)
        fig_V.savefig('plot_results/sim_path.pdf',format='pdf')

    fig_V, ax_V = plt.subplots()
    plt.tight_layout
    curvature_data_sim = curvature_data_sim[np.abs(curvature_data_sim.curvature-curvature_data_sim.curvature.mean())<= (3.0 * curvature_data_sim.curvature.std())]
    x = curvature_data_sim['velocity']
    y = curvature_data_sim['curvature']
    coefficients = np.polyfit(x,y,5)
    poly = np.poly1d(coefficients)
    new_x = np.linspace(min(x), max(x))
    new_y = poly(new_x)
    ax_V.plot(x,y,"o",new_x,new_y)
    ax_V.scatter(curvature_data_sim['velocity'],curvature_data_sim['curvature'],s =0.3)
    ax_V.plot(x,y,"o",new_x,new_y)
    ax_V.set_xlabel("velocity [µm/s]",fontdict=None,labelpad=0, rotation=0,)
    ax_V.set_ylabel("curvature",fontdict=None,labelpad=0)
    ax_V.ticklabel_format(axis='both', style = 'sci')
    ax_V.legend(loc='best',fontsize=8)
    ax_V.grid(True,linestyle='-')
    ax_V.set_title("Simulated Curvature vs. Velocity", fontdict= None, loc = 'center',pad=5)
    fig_V.savefig('plot_results/sim_curv time %i.pdf' %k,format='pdf')

    fig_V, ax_V = plt.subplots()
    plt.tight_layout
    num_bins = 100
    n, bins, patches = ax_V.hist(total_velocities["Velocity"], num_bins, density=1)
    sigma = total_velocities["Velocity"].std()
    mu = total_velocities["Velocity"].mean()
    print("mu")
    print(mu)
    print("sigma")
    print(sigma)
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    ax_V.plot(bins, y, '--')
    ax_V.set_xlabel("velocity [µm/s]",fontdict=None,labelpad=0, rotation=0,fontsize=10)
    ax_V.set_ylabel("Probability Density",fontdict=None,labelpad=0,fontsize=10)
    ax_V.ticklabel_format(axis='both', style = 'sci')
    #ax_V.legend(loc='best',fontsize=8)
    ax_V.grid(True,linestyle='-')
    #ax_V.set_title("Experimental Curvature vs. Velocity", fontdict= None, loc = 'center',pad=5)
    fig_V.tight_layout()
    plt.show()
    fig_V.savefig('plot_results/simulated_total_velocity_histogram time %i_.pdf' %k,format='pdf')
##End : plot curvature of simulated dumbbell data ##




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




fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(position['Position in X']* pow(10,6), position['Position in Y']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.plot(position['Position in X 2'] * pow(10,6), position['Position in Y 2']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
ax_V.plot(position['Position in X 3']* pow(10,6), position['Position in Y 3']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)
ax_V.plot(position['Position in X 4']* pow(10,6), position['Position in Y 4']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)


ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
#ax_V.set_title("Dumbbell Trajectories ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('path.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(rotator_position['Position in X']* pow(10,6), rotator_position['Position in Y']* pow(10,6), marker=next(marker),linestyle=next(line),markersize = 0.8)
ax_V.plot(rotator_position['Position in X 2'] * pow(10,6), rotator_position['Position in Y 2']* pow(10,6), marker=next(marker),linestyle=next(line), markersize = 0.8)

ax_V.set_xlabel("x(t) [µm]",fontdict=None,labelpad=0)
ax_V.set_ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
ax_V.ticklabel_format(axis='both', style = 'sci')
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Rotator Trajectories ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('path.pdf',format='pdf')

cm = plt.cm.get_cmap('RdYlBu')
fig_V, ax_V = plt.subplots()
plt.tight_layout

cm = plt.cm.get_cmap('RdYlBu')
fig_V, ax_V = plt.subplots()
plt.tight_layout

sc = plt.scatter(x=position['Position in X']* pow(10,6),y=position['Position in Y']* pow(10,6),c=position['Time Step'] * 2 * pow(10,-3), cmap=cm, s=0.7)
plt.colorbar(sc)

plt.xlabel("x(t) [µm]",fontdict=None,labelpad=0,fontsize=10)
plt.ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,fontsize=10)
plt.legend(loc='best',fontsize=8)
plt.title("path ", fontdict= None, loc = 'center',pad=5)
plt.savefig('path_with_colors.pdf',format='pdf')


        #ax_V.scatter(position['Position in X'], position['Position in Y'], c = position['Time Step'])

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
ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Rotator Trajectories", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('rot_path.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(position['Time Step'], position['Position in X'], marker=next(marker),linestyle=next(line))
ax_V.plot(position['Time Step'], position['Position in Y'], marker=next(marker),linestyle=next(line))

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
sc = plt.scatter(x=triangle_position['COM Position in X 4']* pow(10,6),y=triangle_position['COM Position in Y 4']* pow(10,6),c=triangle_position['Time Step'] * 2.0 * pow(10,-3), cmap=cm,s = 0.8)
colorb = plt.colorbar(sc)
colorb.ax.get_yaxis().labelpad = 15
colorb.ax.set_ylabel('Time [s]', rotation=270)
plt.xlabel("x(t) [µm]",fontdict=None,labelpad=0)
plt.ylabel("y(t) [µm]",fontdict=None,labelpad=0, rotation=90,)
plt.ticklabel_format(axis='both', style = 'sci')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
ax_V.set_title("Triangle Trajectories", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('tri_color_path.pdf',format='pdf')



fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot(FMA['time step'], FMA['F-M*a'], marker=next(marker),linestyle=next(line))


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

ax_V.plot((v_w_data['Frequency']), v_w_data[' Velocity'] * pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.set_xscale('log')
ax_V.axhline(linewidth=1, color='k', linestyle = '-.')
ax_V.set_xlabel("Frequency [Hz]",fontdict=None,labelpad=0)
ax_V.set_ylabel("Velocity [m/s]",fontdict=None,labelpad=0, rotation=90,)
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
#ax_V.set_title("Velocity vs Frequency in micrometer/s ", fontdict= None, loc = 'center',pad=5)
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

v_w_2 = pd.read_csv("../../v_w_2.csv")
print("v_w_")
print(v_w_2)
fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot((v_w_2['frequency']), v_w_2['U'] * pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.plot((v_w_2['frequency']), v_w_2['U A']* pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.plot((v_w_2['frequency']), v_w_2['U B']* pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.set_xscale('log')
ax_V.axhline(linewidth=1, color='k', linestyle = '-.')
ax_V.set_xlabel("Frequency [Hz]",fontdict=None,labelpad=0, fontsize= 10)
ax_V.set_ylabel("Velocity [m/s]",fontdict=None,labelpad=0, rotation=90,fontsize=10)
ax_V.legend(loc='best',fontsize=16)
L =plt.legend()
L.get_texts()[1].set_text('$U_B$')
L.get_texts()[2].set_text('$U_A$')
ax_V.grid(True,linestyle='-')
#ax_V.set_title("U vs Frequency in micrometer/s ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('Velocity_frequency.pdf',format='pdf')

fig_V, ax_V = plt.subplots()
plt.tight_layout



ax_V.plot((v_w_2['frequency']), v_w_2['U'] * pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.plot((v_w_data['Frequency']), v_w_data[' Velocity'] * pow(10,6), marker=next(marker),linestyle=next(line))
ax_V.set_xscale('log')
ax_V.axhline(linewidth=1, color='k', linestyle = '-.')
ax_V.set_xlabel("Frequency [Hz]",fontdict=None,labelpad=0, fontsize= 10)
ax_V.set_ylabel("Velocity [m/s]",fontdict=None,labelpad=0, rotation=90,fontsize=10)
ax_V.legend(loc='best',fontsize=16)
L =plt.legend()
L.get_texts()[0].set_text('U')
L.get_texts()[1].set_text('V')
ax_V.grid(True,linestyle='-')
#ax_V.set_title("U vs Frequency in micrometer/s ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('Velocity_frequency.pdf',format='pdf')


fig_V, ax_V = plt.subplots()
plt.tight_layout

ax_V.plot((v_w_2['frequency']), v_w_2['RE CM'] , marker=next(marker),linestyle=next(line))
ax_V.plot((v_w_2['frequency']), v_w_2['IM CM']  , marker=next(marker),linestyle=next(line))
ax_V.set_xscale('log')
ax_V.axhline(linewidth=1, color='k', linestyle = '-.')
ax_V.set_xlabel("Frequency [Hz]",fontdict=None,labelpad=0)
ax_V.set_ylabel("Clausius Mossotti factor",fontdict=None,labelpad=0, rotation=90,)
L =plt.legend()
L.get_texts()[0].set_text('Re[CM]')
L.get_texts()[1].set_text('Im[CM]')
#ax_V.legend(loc='best',fontsize=8)
ax_V.grid(True,linestyle='-')
#ax_V.set_title("CM vs Frequency in micrometer/s ", fontdict= None, loc = 'center',pad=5)
fig_V.savefig('CM_frequency.pdf',format='pdf')



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
