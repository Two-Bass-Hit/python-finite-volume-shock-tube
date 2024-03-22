
import numpy as np
import scipy as sp

def WtoFW(W, FW, geom_array, cv, cp, r, i , n): #i = node, n = time (0 = t or 1 = t + 1), 

  F = geom_array[1, i] # node area 

  rho = W[0, i, n] / F
  u = W[1, i, n] / W[0, i, n]
  e0 = W[2, i, n] / W[0, i, n]
  e = e0 - u**2 / 2
  T = e / cv
  P = rho * r * T
  h0 = cp * T + u**2 / 2
  FW[0, i, n] = W[1, i, n]
  FW[1, i, n] = W[1, i, n] * u + P * F
  FW[2, i, n] = h0 * W[1, i, n]


def WtoFW_vector_plus_one_over_2(W, FW, geom_array, cv, cp, r, i): #i = node, 

  F = (geom_array[1, i] + geom_array[1, i + 1]) / 2 # area at i + 1/2

  rho = W[0, 0, 0] / F
  u = W[1, 0, 0] / W[0, 0, 0]
  e0 = W[2, 0, 0] / W[0, 0, 0]
  e = e0 - u**2 / 2
  T = e / cv
  P = rho * r * T
  h0 = cp * T + u**2 / 2
  FW[0, 0, 0] = W[1, 0, 0]
  FW[1, 0, 0] = W[1, 0, 0] * u + P * F
  FW[2, 0, 0] = h0 * W[1, 0, 0]

def WtoFW_vector_minus_one_over_2(W, FW, geom_array, cv, cp, r, i): #i = node, 

  F = (geom_array[1, i] + geom_array[1, i - 1]) / 2 # area at i + 1/2 

  rho = W[0, 0, 0] / F
  u = W[1, 0, 0] / W[0, 0, 0]
  e0 = W[2, 0, 0] / (rho * F)
  e = e0 - u**2 / 2
  T = e / cv
  P = rho * r * T
  h0 = cp * T + u**2 / 2
  FW[0, 0, 0] = W[1, 0, 0]
  FW[1, 0, 0] = W[1, 0, 0] * u + P * F
  FW[2, 0, 0] = h0 * W[1, 0, 0]  