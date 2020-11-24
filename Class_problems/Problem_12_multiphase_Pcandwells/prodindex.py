"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Update wells
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

import numpy as np
from rel_perm import rel_perm

def prodindex(i,well,reservoir,fluid,petro,numerical,Sw):  
    
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
        krw,kro = rel_perm(petro,Sw[kblock,0])
        
        Jwell_w = 6.33E-3*(2*np.pi*perm*krw*reservoir.h)/(fluid.muw[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
        Jwell_o = 6.33E-3*(2*np.pi*perm*kro*reservoir.h)/(fluid.muo[kblock,0]*fluid.Bo[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
     
    elif well.direction[i][0] == 'hx': #horizontal x-axis
        #calculating equivalent radius
        kzky  = reservoir.permz[kblock,0] / reservoir.permy[kblock,0] 
        kykz  = 1.0/kzky
        req   = (0.28*np.sqrt(np.sqrt(kzky)*numerical.dy[kblock,0]**2.0 + np.sqrt(kykz)*reservoir.h**2.0))/ (kzky**0.25 + kykz**0.25)
        #req   = 0.2*numerical.dx[kblock,0] #Equivalent radius Peaceman correction [feet]
        
        #calculating the productivity index
        perm  = np.sqrt(reservoir.permy[kblock,0]*reservoir.permz[kblock,0]) #effective permeability
        krw,kro = rel_perm(petro,Sw[kblock,0])
        
        Jwell_w = 6.33E-3*(2*np.pi*perm*krw*numerical.dx[kblock,0])/(fluid.muw[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
        Jwell_o = 6.33E-3*(2*np.pi*perm*kro*numerical.dx[kblock,0])/(fluid.muo[kblock,0]*fluid.Bo[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))

    elif well.direction[i][0] == 'hy': #horizontal y-axis
        #calculating equivalent radius
        kzkx  = reservoir.permz[kblock,0] / reservoir.permx[kblock,0] 
        kxkz  = 1.0/kzkx
        req   = (0.28*np.sqrt(np.sqrt(kzkx)*numerical.dx[kblock,0]**2.0 + np.sqrt(kxkz)*reservoir.h**2.0))/ (kzkx**0.25 + kxkz**0.25)
        #req   = 0.2*numerical.dx[kblock,0] #Equivalent radius Peaceman correction [feet]
        
        #calculating the productivity index
        perm  = np.sqrt(reservoir.permx[kblock,0]*reservoir.permz[kblock,0]) #effective permeability
        krw,kro = rel_perm(petro,Sw[kblock,0])
        
        Jwell_w = 6.33E-3*(2*np.pi*perm*krw*numerical.dy[kblock,0])/(fluid.muw[kblock,0]*fluid.Bw[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
        Jwell_o = 6.33E-3*(2*np.pi*perm*kro*numerical.dy[kblock,0])/(fluid.muo[kblock,0]*fluid.Bo[kblock,0]*(np.log((req/well.rw[i][0])) + well.skin[i][0]))
        
    return Jwell_w,Jwell_o;