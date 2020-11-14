from scipy.sparse import lil_matrix, csr_matrix

#fluid, reservoir and simulation parameters   
def myarrays(fluid,reservoir,numerical,BC,IC):
    #Setting up matrix T, B, and Q
    T = lil_matrix((numerical.N, numerical.N))
    B = lil_matrix((numerical.N, numerical.N))
    Q = lil_matrix((numerical.N, 1))
    G = lil_matrix((numerical.N, 1))

    trans = reservoir.k[0] * reservoir.W * reservoir.h /(fluid.mu * fluid.Bw * numerical.dx[0]) #Transmissibility vector 
    
    for l in range(0,numerical.N):
        if (l+1) % numerical.Nx != 1:  #not on left boundary     
             T[l,l-1] = -trans                     
             T[l,l]   = T[l,l] - T[l,l-1]       
        else:                     #left boundary  
            if 'Neumann' in BC.type[0]:
                None
            elif 'Dirichlet' in BC.type[0]:
                T[l,l] = T[l,l] + 2*trans
                Q[l,0] = Q[l,0] + 2 * trans * (BC.value[0][0] - (fluid.rho/144.0)*BC.depth[0][0])* 6.33E-3      
        
        if (l+1) % numerical.Nx != 0:  
             T[l,l+1] = -trans                   
             T[l,l]   = T[l,l] - T[l,l+1]       
        else:                     #right  boundary
            if 'Neumann' in BC.type[1]:
                None
            elif 'Dirichlet' in BC.type[1]:
                T[l,l] = T[l,l] + 2 * trans
                Q[l,0] = Q[l,0] + 2 * trans * (BC.value[1][0] - (fluid.rho/144.0)*BC.depth[1][0])* 6.33E-3                 
            
        if int(l / numerical.Nx) > 0:  #not bottom boundary
             T[l,l-numerical.Nx] = - trans
             T[l,l]   = T[l,l] - T[l,l-numerical.Nx] 
        else:                   #bottom boundary
            if 'Neumann' in BC.type[2]:
                None
            elif 'Dirichlet' in BC.type[2]:
                T[l,l] = T[l,l] + 2 * trans
                Q[l,0] = Q[l,0] + 2 * trans * (BC.value[2][0] - (fluid.rho/144.0)*numerical.D[l,0] ) * 6.33E-3      
        
        if int(l / numerical.Nx) < numerical.Ny - 1:  #not top boundary
             T[l,l+numerical.Nx] = - trans
             T[l,l]   = T[l,l] - T[l,l+numerical.Nx] 
        else:                                         #top boundary
            if 'Neumann' in BC.type[3]:
                None
            elif 'Dirichlet' in BC.type[3]:
                T[l,l] = T[l,l] + 2*trans
                Q[l,0] = Q[l,0] + 2 * trans * (BC.value[3][0] - (fluid.rho/144.0)*numerical.D[l,0] ) * 6.33E-3                
        
        B[l,l] = numerical.dx[0,0] * numerical.dy[0,0] * reservoir.h * reservoir.phi[l,0] * fluid.ct / fluid.Bw #accumulation

    T = (6.33E-3 * T).tocsr() #multiplying with the conversion factor 
    G = csr_matrix((fluid.rho/144.0)*(T@numerical.D))
    B = B.tocsr()
    Q = Q.tocsr()
    return T, B, Q, G;
