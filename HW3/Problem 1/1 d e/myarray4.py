"""
reservoir simulation assignment 4
2D reservoir simulation: Arrays
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 11/10/2020
"""

from scipy.sparse import lil_matrix, csr_matrix
import numpy as np
from courant import courant

class numerical:
    def __init__(self):
        self.N  = []
class IC:
    def __init__(self):
        self.P  = []
class fluid:
    def __init__(self):
        self.N  = []
class reservoir:
    def __init__(self):
        self.P  = []

#fluid, reservoir and simulation parameters   
#def inputfile(reservoir,numerical,BC,IC):
def myarrays4(reservoir,fluid,numerical,IC):

    #Setting up matrix T, B, and Q
    ETA = lil_matrix((numerical.N, numerical.N))
    CR  = lil_matrix((numerical.N, numerical.N))
    
    for i in range(0,numerical.N):
        if (i+1) % numerical.Nx != 0:  #not on right boundary
            ETA[i,i+1] = -1*numerical.eta
            ETA[i,i]   =  ETA[i,i] + 1.0*numerical.eta             
            
            if IC.P[i,0] > IC.P[i+1,0]:
                CR[i,i]  = CR[i,i] + 1.0*courant(i,i+1,IC.P,numerical,fluid,reservoir)
            elif IC.P[i,0] <= IC.P[i+1,0]:
                CR[i,i+1]= CR[i,i+1] - 1.0*courant(i,i+1,IC.P,numerical,fluid,reservoir)
                   
        if (i+1) % numerical.Nx != 1:  #not on left boundary  
            ETA[i,i-1] = -1*numerical.eta
            ETA[i,i]   =  ETA[i,i] + 1.0*numerical.eta    

            if IC.P[i,0] < IC.P[i-1,0]:
                CR[i,i-1]  = -1.0*courant(i,i-1,IC.P,numerical,fluid,reservoir)
            elif IC.P[i,0] >= IC.P[i-1,0]:
                CR[i,i]    =   CR[i,i] + 1.0*courant(i,i-1,IC.P,numerical,fluid,reservoir)
       
        if int(i / numerical.Nx) > 0:  #not bottom boundary
            ETA[i,i-numerical.Nx] = -1.0*numerical.eta 
            ETA[i,i]   =  ETA[i,i] + 1.0*numerical.eta   

            if IC.P[i,0] < IC.P[i-numerical.Nx,0]:
                CR[i,i-numerical.Nx]  = -1.0*courant(i,i-numerical.Nx,IC.P,numerical,fluid,reservoir)
            elif IC.P[i,0] >= IC.P[i-numerical.Nx,0]:
                CR[i,i]    =   CR[i,i] + 1.0*courant(i,i-numerical.Nx,IC.P,numerical,fluid,reservoir)
                   
        if int(i / numerical.Nx) < numerical.Ny - 1:  #not top boundary            
            ETA[i,i+numerical.Nx] = -1.0*numerical.eta 
            ETA[i,i]   =  ETA[i,i] + 1.0*numerical.eta   

            if IC.P[i,0] < IC.P[i+numerical.Nx,0]:
                CR[i,i+numerical.Nx]  = -1.0*courant(i,i+numerical.Nx,IC.P,numerical,fluid,reservoir)
            elif IC.P[i,0] >= IC.P[i+numerical.Nx,0]:
                CR[i,i]    =   CR[i,i] + 1.0*courant(i,i+numerical.Nx,IC.P,numerical,fluid,reservoir)
             
    return ETA,CR;