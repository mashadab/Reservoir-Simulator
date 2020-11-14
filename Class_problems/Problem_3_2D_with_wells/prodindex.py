import numpy as np

#appending the update the wells
def prodindex(i,well,reservoir,fluid,numerical):  
    # using the determined the block wells 
    kblock = well.block[i,0]

    #calculating the productivity index
    req   = 0.2*numerical.dx[kblock,0] #Equivalent radius [feet]
    perm  = np.sqrt(reservoir.permx[kblock,0]*reservoir.permy[kblock,0]) #effective permeability
    Jwell = 6.33E-3*(2*np.pi*perm*reservoir.h)/(fluid.mu[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
    return Jwell;
