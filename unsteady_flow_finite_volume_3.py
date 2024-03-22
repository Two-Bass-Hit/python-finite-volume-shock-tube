# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:35:47 2019

@author: jerome.karlovsky
"""
#can deal with area change in pipes

import numpy as np
import scipy as sp
import func_WtoFW_2
import func_WtoCW
import func_calc_timestep
import func_calc_W_plus_1
import matplotlib.pyplot as plt
import math as math

pipe_length = 1
num_nodes = 101
diameter_left = 0.1
diameter_right = 0.1
dd_by_dx = (diameter_right - diameter_left) / pipe_length
W = np.empty(shape = (5, num_nodes, 2)) #makes empty 3d array (random number not zeros), if printed this is 5 arrays of 101 rows and 2 columns
FW = np.empty(shape = (3, num_nodes, 2)) #makes empty 3d array (random number not zeros), if printed this is 3 arrays of 101 rows and 2 columns
CW = np.empty(shape = (3, num_nodes, 2)) #makes empty 3d array (random number not zeros), if printed this is 3 arrays of 101 rows and 2 columns
geom_array =np.empty(shape = (3, num_nodes)) 

for index in range(0, num_nodes):

  geom_array[0, index] = diameter_left + index / (num_nodes - 1) * (diameter_right - diameter_left) #diameter of node
  geom_array[1, index] = math.pi / 4 * geom_array[0, index] ** 2                              #area of node
  geom_array[2, index] = math.pi / 2 * geom_array[0, index] * dd_by_dx                        # dA/dx at node

dx = pipe_length / (num_nodes - 1) #mesh length

r = 287
cp = 1005
cv = 718
g = cp / cv
#left hand initial pipe properties

P = 500000
T = 1200
rho = P / (r * T)
u = 0
e0 = cv * T
i = 0

while i < 51:

    F = geom_array[1, i]   

    W[0, i, 0] = rho * F      #sets continuity term in nodes 0 - 49
    W[1, i, 0] = rho * u * F  #sets momentum term in nodes 0 - 49
    W[2, i, 0] = rho * e0 * F #sets energy term in nodes 0 - 49
    W[3, i, 0] = P            #sets pressure in nodes 0 - 49
    W[4, i, 0] = T            #sets pressure in nodes 0 - 49
    
    i = i + 1

#right hand initial pipe properties

P = 100000.000
T = 300.000
rho = P / (r * T)
u = 0
e0 = cv * T

#i = 0 #this sets all the pipe nodes to have the same property

while i < 101:

    F = geom_array[1, i] 

    W[0, i, 0] = rho * F      #sets continuity term in the left hand side of pipe in nodes 50 - 101
    W[1, i, 0] = rho * u * F  #sets momentum term in the left hand side of pipe in nodes 50 - 101
    W[2, i, 0] = rho * e0 * F #sets energy term in the left hand side of pipe in nodes 50 - 101
    W[3, i, 0] = P            #sets pressure in the left hand side of pipe in nodes 50 - 101
    W[4, i, 0] = T            #sets temperature in the left hand side of pipe in nodes 50 - 101
    
    i = i + 1

print(W)

i = 0

while i < 101: #initial conditions. This isn't really needed, FW could be calculated in the W-Plus_1 function, but it is nice to have all values available from the main program if need be.
     
    func_WtoFW_2.WtoFW(W, FW, geom_array, cv, cp, r, i, 0)
    func_WtoCW.WtoCW(W, CW, geom_array, i, 0)

    i = i + 1

print(FW)    

time = 0

while time <= 5 * 10**-4:

    dt = 10 #far bigger than it will ever be calculated to be
    #do for each pipe although there is only one pipe for the shock tube test
    dt = func_calc_timestep.calc_timestep(W, dx, dt, g, r) #calcs minimum time for wave 

    courant_number = 0.8
    dt = courant_number * dt
    print(dt)    

    i = 1 #start from 2nd node in

    while i < 100: #loops through all the nodes except for first and last

        func_calc_W_plus_1.calc_W_plus_1(i, W, FW, CW, dt, dx, geom_array, cv, cp, r, courant_number, num_nodes) #includes flux limiter
        
        poo = W[0, i, 1]
        wee = FW[0, i, 1]

        i = i + 1
 
    #moves properties from time t + 1 to time t ready for use in the next timestep
    W[:, 1:100, 0] = W[:,1:100, 1] #moves data in nodes 1 to 99 at t + 1 back to t for next timestep
    FW[:, 1:100, 0] = FW[:,1:100, 1] #moves data in nodes 1 to 99 at t + 1 back to t for next timestep
    CW[:, 1:100, 0] = CW[:,1:100, 1] #moves data in nodes 1 to 99 at t + 1 back to t for next timestep

    time = time + dt

P_array = W[3, :, 0]
T_array = W[4, :, 0]   

print(P_array)
print(T_array)

fig = plt.figure(figsize=(12,3))
ax1 = fig.add_subplot(111)
ax1.plot(P_array,  linewidth=.5)

fig = plt.figure(figsize=(12,3))
ax2 = fig.add_subplot(111)
ax2.plot(T_array,  linewidth=.5)

plt.show() 