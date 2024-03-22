
import numpy as np
import scipy as sp



def WtoFW(W, FW, F, cv, cp, r):

  rho = W[0] / F
  u = W[1] / W[0]
  e0 = W[2] / (rho * F)
  e = e0 - u**2 / 2
  T = e / cv
  P = rho * r * T
  h0 = cp * T + u**2 / 2
  FW[0] = W[1]
  FW[1] = W[1] * u + P * F
  FW[2] = h0 * W[1]

