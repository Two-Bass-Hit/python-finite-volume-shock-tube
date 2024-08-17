import numpy as np
import scipy as sp

def flux_limiter_tvd(i, W, courant_number, num_nodes, geom_array): 

  delta_W_i_plus_1_over_2 = np.empty(shape = (3)) #creates 3x1 vector
  delta_W_i_minus_1_over_2 = np.empty(shape = (3))
  delta_W_i_plus_3_over_2 = np.empty(shape = (3))
  delta_W_i_minus_3_over_2 = np.empty(shape = (3))

#calculates delta Wn 

  delta_W_i_plus_1_over_2[:3] = W[:3, i + 1, 0] / geom_array[1, i + 1] - W[:3, i, 0] / geom_array[1, i]
  delta_W_i_minus_1_over_2[:3] =W[:3, i, 0] / geom_array[1, i] - W[:3, i - 1, 0] / geom_array[1, i - 1]

  if i == 1:

    delta_W_i_plus_3_over_2[:3] = W[:3, i + 2, 0] / geom_array[1, i + 2] - W[:3, i + 1, 0] / geom_array[1, i + 1]
    delta_W_i_minus_3_over_2[:3] = delta_W_i_minus_1_over_2[:3] #assume gradient from virtual node 0-1 to node 0 is the same as from node 0 to 1

  elif i == num_nodes - 2:  #last node that will be calculated will be second last node

    delta_W_i_minus_3_over_2[:3] = W[:3, i - 1, 0] / geom_array[1, i - 1] - W[:3, i - 2, 0] / geom_array[1, i - 2]
    delta_W_i_plus_3_over_2[:3] = delta_W_i_plus_1_over_2[:3]#assume gradient from last node to virtual node last+1 is the same as from node last-1 to last
  
  else:
    
    delta_W_i_plus_3_over_2[:3] = W[:3, i + 2, 0] / geom_array[1, i + 2] - W[:3, i + 1, 0] / geom_array[1, i + 1]
    delta_W_i_minus_3_over_2[:3] = W[:3, i - 1, 0] / geom_array[1, i - 1] - W[:3, i - 2, 0] / geom_array[1, i - 2] 
  
#calculates r values

  r_i_minus_1_plus = np.inner(delta_W_i_minus_3_over_2, delta_W_i_minus_1_over_2) / np.inner(delta_W_i_minus_1_over_2, delta_W_i_minus_1_over_2)
  r_i_plus = np.inner(delta_W_i_minus_1_over_2, delta_W_i_plus_1_over_2) / np.inner(delta_W_i_plus_1_over_2, delta_W_i_plus_1_over_2)
  r_i_minus = np.inner(delta_W_i_minus_1_over_2, delta_W_i_plus_1_over_2)/np.inner(delta_W_i_minus_1_over_2, delta_W_i_minus_1_over_2)
  r_i_plus_1_minus = np.inner(delta_W_i_plus_1_over_2, delta_W_i_plus_3_over_2)/np.inner(delta_W_i_plus_1_over_2, delta_W_i_plus_1_over_2)

  #calculates C

  if courant_number <= 0.5:

    C = courant_number * (1 - courant_number)

  else:
    
    C = 0.25   

  #calculates flux limiter
  
  G1 = 0.5 * C * (1 - calc_phi(r_i_plus))
  G2 = 0.5 * C * (1 - calc_phi(r_i_plus_1_minus))
  G3 = 0.5 * C * (1 - calc_phi(r_i_minus_1_plus))
  G4 = 0.5 * C * (1 - calc_phi(r_i_minus))
  
  TVD_scheme = (G1 + G2) * delta_W_i_plus_1_over_2 * (geom_array[1, i] + geom_array[1, i + 1]) / 2 - (G3 + G4) * delta_W_i_minus_1_over_2 * (geom_array[1, i] + geom_array[1, i - 1]) / 2

  return TVD_scheme  





def calc_phi(r):

  if r > 0:

    phi = min(2 * r, 1)  

  else:

    phi = 0  

  return phi
