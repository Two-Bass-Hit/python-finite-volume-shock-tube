
import numpy as np
import scipy as sp
import func_WtoFW_2
import func_flux_limiter_tvd
import func_WtoCW
import math as math
import func_W_add_T_and_P

def calc_W_plus_1(i, W, FW, CW, dt, dx, geom_array, cv, cp, r, courant_number, num_nodes): #i = node in base zero

    W_plus_half = np.empty(shape = (5, 1, 1))
    W_minus_half = np.empty(shape = (5, 1, 1))
    FW_plus_half = np.empty(shape = (3, 1, 1))
    FW_minus_half = np.empty(shape = (3, 1, 1))
    CW_plus_half = np.empty(shape = (3, 1, 1))
    CW_minus_half = np.empty(shape = (3, 1, 1))

    #First step lax Friedrich method loops through the 3 terms in vector
    
    

    W_plus_half[:3, 0, 0] = 0.5 * (W[:3, i + 1, 0] + W[:3, i, 0]) - dt / (2 * dx) * (FW[:3, i + 1, 0] - FW[:3, i, 0]) - dt / 4 * (CW[:3, i + 1, 0] + CW[:3, i, 0])
    W_minus_half[:3, 0, 0] = 0.5 * (W[:3, i, 0] + W[:3, i - 1, 0]) - dt / (2 * dx) * (FW[:3, i, 0] - FW[:3, i - 1, 0]) - dt / 4 * (CW[:3, i, 0] + CW[:3, i - 1, 0])
    
    F_plus_half = (geom_array[1, i] + geom_array[1, i + 1]) / 2  # node area at i + 1/2 calculated as average of the end areas not from diameter
    F_minus_half = (geom_array[1, i] + geom_array[1, i - 1]) / 2   # node area at i - 1/2 calculated as average of the end areas not from diameter
   
    func_W_add_T_and_P.W_add_T_and_P(W_plus_half, F_plus_half, cv, cp, r, 0, 0)  #calculates T and P at time t + 1/2 and position i + 1/2 and adds to W_plus_half
    func_W_add_T_and_P.W_add_T_and_P(W_minus_half, F_minus_half, cv, cp, r, 0, 0) #calculates T and P at time t - 1/2 and position i - 1/2 and adds to W_minus_half

    func_WtoFW_2.WtoFW_vector_plus_one_over_2(W_plus_half, FW_plus_half, geom_array, cv, cp, r, i)   #calculates FW plus half
    func_WtoFW_2.WtoFW_vector_minus_one_over_2(W_minus_half, FW_minus_half, geom_array, cv, cp, r, i) #calculates FW minus half
    
    func_WtoCW.WtoCW_vector_plus_1_over_2(W_plus_half, CW_plus_half, geom_array, i)   #calculates CW plus half
    func_WtoCW.WtoCW_vector_minus_1_over_2(W_minus_half, CW_minus_half, geom_array, i) #calculates CW minus half


    #second step leapfrog method loops through the 3 terms and fills in the next time position in the state array

    W[:3, i, 1] = W[:3, i, 0] - dt / dx * (FW_plus_half[:3, 0, 0] - FW_minus_half[:3, 0, 0]) - dt / 2 * (CW_plus_half[:3, 0, 0] + CW_minus_half[:3, 0, 0])
    
    flux_lim = func_flux_limiter_tvd.flux_limiter_tvd(i, W, courant_number, num_nodes, geom_array)

    W[:3, i, 1] = W[:3, i, 1] + flux_lim[:3]

    F = geom_array[1, i] # node area 
    rho = W[0, i, 1] / F
    u = W[1, i, 1] / W[0, i, 1]
    e0 = W[2, i, 1] / (rho * F)
    e = e0 - u**2 / 2
    T = e / cv
    P = rho * r * T   

    W[3, i, 1] = P
    W[4, i, 1] = T

    func_WtoFW_2.WtoFW(W, FW, geom_array, cv, cp, r, i, 1) # fills in the flux aray at time t + 1 using the state array at t + 1, i = node, 1 denotes time t + 1 in state and flux array
    func_WtoCW.WtoCW(W, CW, geom_array, i, 1) # fills in the cource aray at time t + 1 using the state array at t + 1, i = node, 1 denotes time t + 1 in state and source array