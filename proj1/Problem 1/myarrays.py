"""
reservoir simulation project 1 (2020): Problem 1
2D reservoir simulation: Building arrays
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

from scipy.sparse import lil_matrix, csr_matrix
from Thalf import Thalf   #for calculating transmissibility

#fluid, reservoir and simulation parameters   
def myarrays(fluid,reservoir,numerical,BC,IC):
    #Setting up matrix T, B, and Q
    T = lil_matrix((numerical.N, numerical.N))
    B = lil_matrix((numerical.N, numerical.N))
    Q = lil_matrix((numerical.N, 1))
    G = lil_matrix((numerical.N, 1))

    #trans = reservoir.k[0,0] * reservoir.W * reservoir.h /(fluid.mu[0,0] * fluid.Bw[0,0] * numerical.dx[0,0]) #Transmissibility vector 
    
    for l in range(0,numerical.N):
        if (l+1) % numerical.Nx != 1:  #not on left boundary     
             T[l,l-1] = -Thalf(l,l-1,'x',fluid,reservoir,numerical)                   
             T[l,l]   = T[l,l] - T[l,l-1]       
        else:                     #left boundary  
            if 'Neumann' in BC.type[0]:
                None
            elif 'Dirichlet' in BC.type[0]:
                T[l,l] = T[l,l] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)* (BC.value[0][0])* 6.33E-3      #without gravity  
                Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)* (BC.value[0][0] - (fluid.rho/144.0)*numerical.D[l,0] )* 6.33E-3 #with gravity   
        
        if (l+1) % numerical.Nx != 0:  #not on right boundary
             T[l,l+1] = -Thalf(l,l+1,'x',fluid,reservoir,numerical)                                    
             T[l,l]   = T[l,l] - T[l,l+1]       
        else:                     #right  boundary
            if 'Neumann' in BC.type[1]:
                None
            elif 'Dirichlet' in BC.type[1]:
                T[l,l] = T[l,l] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical) * (BC.value[1][0])* 6.33E-3  #without gravity                 
                Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical) * (BC.value[1][0] - (fluid.rho/144.0)*numerical.D[l,0] )* 6.33E-3  #with gravity  
    
        if int(l / numerical.Nx) > 0:  #not bottom boundary
             T[l,l-numerical.Nx] = - Thalf(l,l-numerical.Nx,'y',fluid,reservoir,numerical)
             T[l,l]   = T[l,l] - T[l,l-numerical.Nx] 
        else:                   #bottom boundary
            if 'Neumann' in BC.type[2]:
                None
            elif 'Dirichlet' in BC.type[2]:
                T[l,l] = T[l,l] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[2][0]) * 6.33E-3 #without gravity 
                Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[2][0] - (fluid.rho/144.0)*numerical.D[l,0] ) * 6.33E-3 #with gravity 
        
        if int(l / numerical.Nx) < numerical.Ny - 1:  #not top boundary
             T[l,l+numerical.Nx] = - Thalf(l,l+numerical.Nx,'y',fluid,reservoir,numerical)
             T[l,l]   = T[l,l] - T[l,l+numerical.Nx] 
        else:                                         #top boundary
            if 'Neumann' in BC.type[3]:
                None
            elif 'Dirichlet' in BC.type[3]:
                T[l,l] = T[l,l] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[3][0]) * 6.33E-3 #without gravity                 
                Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[3][0] - (fluid.rho/144.0)*numerical.D[l,0] ) * 6.33E-3 #with gravity   
    
        B[l,l] = numerical.dx[l,0] * numerical.dy[l,0] * reservoir.h * reservoir.phi[l,0] * fluid.ct / fluid.Bw[l,0] #accumulation

    T = (6.33E-3 * T).tocsr() #multiplying with the conversion factor 
    G = csr_matrix((fluid.rho/144.0)*(T@numerical.D))
    B = B.tocsr()
    Q = Q.tocsr()
    return T, B, Q, G;
