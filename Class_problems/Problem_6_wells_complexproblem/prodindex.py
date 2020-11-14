import numpy as np
'''
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
'''

def prodindex(i,well,reservoir,fluid,numerical):  
    
    # using the determined the block wells 
    kblock = well.block[i,0]

    if well.direction[i][0] == 'v':
        #calculating equivalent radius
        kykx  = reservoir.permy[kblock,0] / reservoir.permx[kblock,0] 
        kxky  = 1.0/kykx
        req   = (0.28*np.sqrt(np.sqrt(kykx)*numerical.dx[kblock,0]**2.0 + np.sqrt(kxky)*numerical.dy[kblock,0]**2.0))/ (kykx**0.25 + kxky**0.25)
        #req   = 0.2*numerical.dx[kblock,0] #Equivalent radius Peaceman correction [feet]
        
        #calculating the productivity index
        perm  = np.sqrt(reservoir.permx[kblock,0]*reservoir.permy[kblock,0]) #effective permeability
        Jwell = 6.33E-3*(2*np.pi*perm*reservoir.h)/(fluid.mu[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
    
    elif well.direction[i][0] == 'hx': #horizontal x-axis
        #calculating equivalent radius
        kzky  = reservoir.permz[kblock,0] / reservoir.permy[kblock,0] 
        kykz  = 1.0/kzky
        req   = (0.28*np.sqrt(np.sqrt(kzky)*numerical.dy[kblock,0]**2.0 + np.sqrt(kykz)*reservoir.h**2.0))/ (kzky**0.25 + kykz**0.25)
        #req   = 0.2*numerical.dx[kblock,0] #Equivalent radius Peaceman correction [feet]
        
        #calculating the productivity index
        perm  = np.sqrt(reservoir.permy[kblock,0]*reservoir.permz[kblock,0]) #effective permeability
        Jwell = 6.33E-3*(2*np.pi*perm*numerical.dx[kblock,0])/(fluid.mu[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
              
    return Jwell;