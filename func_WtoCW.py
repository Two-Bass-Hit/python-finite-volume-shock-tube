import numpy as np
import scipy as sp

def WtoCW(W, CW, geom_array, i, n): #i = node, n = time (0 = t or 1 = t + 1), 
  
  CW[0, i, n] = 0
  CW[1, i, n] = - W[3, i, n] * geom_array[2, i]  #-p * dA/dx
  CW[2, i, n] = 0  



def WtoCW_vector_minus_1_over_2(W, CW, geom_array, i): #i = node

  CW[0, 0, 0] = 0
  CW[1, 0, 0] = - W[3, 0, 0] * (geom_array[2, i] + geom_array[2, i - 1]) / 2 #-p * dA/dx
  CW[2, 0, 0] = 0


def WtoCW_vector_plus_1_over_2(W, CW, geom_array, i): #i = node

  CW[0, 0, 0] = 0
  CW[1, 0, 0] = - W[3, 0, 0] * (geom_array[2, i] + geom_array[2, i + 1]) / 2 #-p * dA/dx
  CW[2, 0, 0] = 0
  