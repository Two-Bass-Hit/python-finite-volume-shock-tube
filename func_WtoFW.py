
import numpy as np
import scipy as sp



def WtoFW(W, FW, F, cv, cp, r,  i , n): #i = node, n = time (0 = t or 1 = t + 1), F = area

  rho = W[0, i, n] / F
  u = W[1, i, n] / W[0, i, n]
  e0 = W[2, i, n] / (rho * F)
  e = e0 - u**2 / 2
  T = e / cv
  P = rho * r * T
  h0 = cp * T + u**2 / 2
  FW[0, i, n] = W[1, i, n]
  FW[1, i, n] = W[1, i, n] * u + P * F
  FW[2, i, n] = h0 * W[1, i, n]

