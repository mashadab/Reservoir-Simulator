"""
reservoir simulation project 1 (2020): Problem 1
2D reservoir simulation: Producer index
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

import numpy as np

#appending the update the wells
def prodindex(i,well,reservoir,fluid,numerical):  
    # using the determined the block wells 
    kblock = well.block[i,0]

    #calculating equivalent radius
    kykx  = reservoir.permy[kblock,0] / reservoir.permx[kblock,0] 
    kxky  = 1.0/kykx
    req   = (0.28*np.sqrt(np.sqrt(kykx)*numerical.dx[kblock,0]**2.0 + np.sqrt(kxky)*numerical.dy[kblock,0]**2.0))/ (kykx**0.25 + kxky**0.25)
    #req   = 0.2*numerical.dx[kblock,0] #Equivalent radius Peaceman correction [feet]
    
    #calculating the productivity index
    perm  = np.sqrt(reservoir.permx[kblock,0]*reservoir.permy[kblock,0]) #effective permeability
    Jwell = 6.33E-3*(2*np.pi*perm*reservoir.h)/(fluid.mu[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
    
    return Jwell;
