from scipy.sparse import lil_matrix, csr_matrix
import numpy as np
from prodindex import prodindex

#appending the update the wells
def updatewells(reservoir,fluid,numerical,P,well,Q): 
    #Setting up J matrix
    J = lil_matrix((numerical.N, numerical.N))
    well.block = np.zeros((len(well.x), 1),dtype='int64')
    well.Jvec  = np.zeros((len(well.x), 1),dtype='float64')    
    for k in range(0,len(well.x)):
        iblock = 0
        for i in range(0,numerical.Nx):
            if well.x[k][0]<(numerical.xc[i,0]+numerical.dx1[i,0]/2) and well.x[k][0]>=(numerical.xc[i,0]-numerical.dx1[i,0]/2)  :
                iblock= i
                break
        jblock = 0
        for j in range(0,numerical.Ny):
            if well.y[k][0]<(numerical.yc[j,0]+numerical.dy1[j,0]/2) and well.y[k][0]>=(numerical.yc[j,0]-numerical.dy1[j,0]/2)  :
                jblock= j
                break
        kblock = iblock + jblock * numerical.Nx
        well.block[k,0] = kblock
        well.Jvec[k,0]  = prodindex(k,well,reservoir,fluid,numerical)
        
        if well.type[k][0] == 2:  #for BHP [psi]
            J[kblock,kblock] = J[kblock,kblock] + well.Jvec[k,0]
            Q[kblock,0]      = Q[kblock,0] + J[kblock,kblock]*well.constraint[k][0]
        elif well.type[k][0] == 1:#for rate [scf/day]  
            Q[kblock,0]      = Q[kblock,0] + well.constraint[k][0] 
    J = J.tocsr()
    return J, Q;
