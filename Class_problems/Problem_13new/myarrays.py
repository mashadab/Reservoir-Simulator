"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Making arrays
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import inv
from Thalf import Thalf   #for calculating transmissibility
from cap_press import cap_press

#fluid, reservoir and simulation parameters   
def myarrays(fluid,reservoir,petro,numerical,BC,P,Pw,Sw,Sw_hyst):
    #Setting up matrix T, B, and Q
    T = lil_matrix((numerical.N, numerical.N))
    Tw= lil_matrix((numerical.N, numerical.N))
    To= lil_matrix((numerical.N, numerical.N))    
    B = lil_matrix((numerical.N, numerical.N))
    d11 = lil_matrix((numerical.N, numerical.N))
    d12 = lil_matrix((numerical.N, numerical.N))  
    d21 = lil_matrix((numerical.N, numerical.N))
    d22 = lil_matrix((numerical.N, numerical.N))  
    D = lil_matrix((numerical.N, numerical.N))
    G = lil_matrix((numerical.N, 1))

    for l in range(0,numerical.N):
        if (l+1) % numerical.Nx != 1:  #not on left boundary     
             Twhalf, Tohalf = Thalf(l,l-1,'x',fluid,reservoir,petro,numerical,P,Sw)
            
             Tw[l,l-1] =-Twhalf                                   
             Tw[l,l]   = Tw[l,l] - Tw[l,l-1]  
            
             To[l,l-1] =-Tohalf                                   
             To[l,l]   = To[l,l] - To[l,l-1]
        else:                     #left boundary  
            if 'Neumann' in BC.type[0]:
                None
            elif 'Dirichlet' in BC.type[0]:
                None
                #T[l,l] = T[l,l] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)* (BC.value[0][0])* 6.33E-3      #without gravity  
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)* (BC.value[0][0] - (fluid.rho/144.0)*numerical.D[l,0] )* 6.33E-3 #with gravity   
        
        if (l+1) % numerical.Nx != 0:  #not on right boundary

             Twhalf, Tohalf = Thalf(l,l+1,'x',fluid,reservoir,petro,numerical,P,Sw)
            
             Tw[l,l+1] =-Twhalf                                   
             Tw[l,l]   = Tw[l,l] - Tw[l,l+1]  
             
             To[l,l+1] =-Tohalf                                   
             To[l,l]   = To[l,l] - To[l,l+1]  
        else:                     #right  boundary
            if 'Neumann' in BC.type[1]:
                None
            elif 'Dirichlet' in BC.type[1]:
                None
                #T[l,l] = T[l,l] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical) * (BC.value[1][0])* 6.33E-3  #without gravity                 
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'x',fluid,reservoir,numerical) * (BC.value[1][0] - (fluid.rho/144.0)*numerical.D[l,0] )* 6.33E-3  #with gravity  
    
        if int(l / numerical.Nx) > 0:  #not bottom boundary
             Twhalf, Tohalf = Thalf(l,l-numerical.Nx,'y',fluid,reservoir,petro,numerical,P,Sw)

             Tw[l,l-numerical.Nx] =-Twhalf                                   
             Tw[l,l]   = Tw[l,l] - Tw[l,l-numerical.Nx] 
             
             To[l,l-numerical.Nx] =-Tohalf                                   
             To[l,l]   = To[l,l] - To[l,l-numerical.Nx]             
             
        else:                   #bottom boundary
            if 'Neumann' in BC.type[2]:
                None
            elif 'Dirichlet' in BC.type[2]:
                None
                #T[l,l] = T[l,l] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[2][0]) * 6.33E-3 #without gravity 
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[2][0] - (fluid.rho/144.0)*numerical.D[l,0] ) * 6.33E-3 #with gravity 

        if int(l / numerical.Nx) < numerical.Ny - 1:  #not top boundary
             Twhalf, Tohalf = Thalf(l,l+numerical.Nx,'y',fluid,reservoir,petro,numerical,P,Sw)

             Tw[l,l+numerical.Nx] =-Twhalf                                   
             Tw[l,l]   = Tw[l,l] - Tw[l,l+numerical.Nx] 
             
             To[l,l+numerical.Nx] =-Tohalf                                   
             To[l,l]   = To[l,l] - To[l,l+numerical.Nx]  
        else:                                         #top boundary
            if 'Neumann' in BC.type[3]:
                None
            elif 'Dirichlet' in BC.type[3]:
                None
                #T[l,l] = T[l,l] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical)
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[3][0]) * 6.33E-3 #without gravity                 
                #Q[l,0] = Q[l,0] + 2 * Thalf(l,l,'y',fluid,reservoir,numerical) * (BC.value[3][0] - (fluid.rho/144.0)*numerical.D[l,0] ) * 6.33E-3 #with gravity   
    
        #B[l,l] = numerical.dx[l,0] * numerical.dy[l,0] * reservoir.h * reservoir.phi[l,0] * fluid.ct / fluid.Bw[l,0] #accumulation
        Vp = numerical.dx[l,0] * numerical.dy[l,0] * reservoir.h * reservoir.phi[l,0]
        Pc,Pcprime = cap_press(petro,Sw[l,0],Sw_hyst[l,0]) 
        d11[l,l] = Vp * Sw[l,0] * (fluid.cw + reservoir.cfr)/(fluid.Bw[l,0] * numerical.dt)
        d12[l,l] = Vp / (fluid.Bw[l,0] * numerical.dt)*(1.0 - Sw[l,0] * reservoir.phi[l,0] * fluid.cw * Pcprime )
        d21[l,l] = Vp * (1-Sw[l,0]) * (fluid.co + reservoir.cfr)/(fluid.Bo[l,0] * numerical.dt)
        d22[l,l] =-Vp /(fluid.Bo[l,0] * numerical.dt)
        
        print(l,Vp,fluid.Bw[l,0],numerical.dt,Sw[l,0],reservoir.phi[l,0],fluid.cw,Pcprime )
        
        D[l,l]   =-(d22[l,l]*d11[l,l]/d12[l,l]) + d21[l,l] 

    d22 = d22.tocsr()
    d12 = d12.tocsr() 
    d21 = d21.tocsr()
    d11 = d11.tocsr() 
    
    Tw= (6.33E-3 * Tw).tocsr()#multiplying with the conversion factor 
    To= (6.33E-3 * To).tocsr()#multiplying with the conversion factor 
    T = (-d22 @ (inv(d12)) @ Tw) + To  #Weighing using the formula given in the sheet                        
    G = -d22 @ inv(d12) @ (Tw @ (P - Pw)) +((-d22 @ inv(d12) @ Tw) * fluid.rhow /144.0 + fluid.rhoo[0,0]/144.0 * To) @ numerical.D
    
    return Tw, To, T, d11, d12, d21, d22, D, G, Pc;
