
def W_add_T_and_P(W, F, cv, cp, r, i, n):  #i = node number, n = time level 0 or 1

    rho = W[0, i, n] / F
    u = W[1, i, n] / W[0, i, n]
    e0 = W[2, i, n] / (rho * F)
    e = e0 - u**2 / 2
    T = e / cv
    P = rho * r * T   

    W[3, i, n] = P
    W[4, i, n] = T