"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Update wells
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

from scipy.sparse import lil_matrix, csr_matrix
import numpy as np
from prodindex import prodindex
from Thalf import Thalf       #for calculating transmissibility
from rel_perm import rel_perm #relative permeability

#appending the update the wells
def updatewells(reservoir,fluid,numerical,petro,P,Sw,well): 
    #Setting up J matrix
    J = lil_matrix((numerical.N, numerical.N))
    well.block = np.zeros((len(well.x), 1),dtype='int64')
    #well.Jvec  = np.zeros((len(well.x), 1),dtype='float64')  
    Qo = lil_matrix((numerical.N, 1))
    Qw = lil_matrix((numerical.N, 1))

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
        #well.Jvec[k,0]  = prodindex(k,well,reservoir,fluid,numerical)
        
        if well.type[k][0] == 2:  #for BHP [psi]
            None
            #J[kblock,kblock] = J[kblock,kblock] + well.Jvec[k,0]
            #Q[kblock,0]      = Q[kblock,0] + J[kblock,kblock]*well.constraint[k][0]
        elif well.type[k][0] == 1:#for rate [scf/day]  
            if well.constraint[k][0] > 0:
                Qw[kblock,0]      = Qw[kblock,0] + well.constraint[k][0] 
                Qo[kblock,0]      = Qo[kblock,0] + 0.0 
            else:
                krw,kro = rel_perm(petro,Sw[kblock,0])
                M       = (kro * fluid.Bw[kblock,0]) / (krw * fluid.Bo[kblock,0])
                fw      = 1.0/(1.0 + M)
                Qw[kblock,0]      = Qw[kblock,0] + fw * well.constraint[k][0] 
                Qo[kblock,0]      = Qo[kblock,0] + (1 - fw) * well.constraint[k][0]                 

    J = J.tocsr()
    Qo= Qo.tocsr()
    Qw= Qw.tocsr()
    return well, Qw,Qo,J;
