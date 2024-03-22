def calc_timestep(W, dx, dt, g, r):
    
    num_nodes = W.shape[1] #length of 2nd dimension of array
    i = 0
    dt_min = dt

    while i <= num_nodes - 1: #remember base 0

        T = W[4, i, 0]
        u = W[1, i, 0] / W[0, i, 0]
        a = (g * r * T)**(1 / 2)
        c = a + abs(u) #propagation velocity
        
        dt_local = dx / c

        if dt_local < dt_min:

            dt_min = dt_local 
      
        i = i + 1
    
    return dt_min