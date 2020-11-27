"""
reservoir simulation assignment 4
1D reservoir simulation Q1.4: Input files
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, identity

class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = [] 
class grid:
    def __init__(self):
        self.xmin = []
class BC:
    def __init__(self):
        self.xmin = []
class IC:
    def __init__(self):
        self.xmin = []
class fluid:
    def __init__(self):
        self.xmin = []

#fluid, reservoir and simulation parameters   
def inputfile():

    #reservoir and fluid properties    
    reservoir.perm   = 1.0     ####################################
    fluid.visc       = 1.0     ####################################    
    
    #Dimensionless properties
    N_Pe             = 1000
    numerical.type   = 'profile'        #Type 'profile' or 'exit'
    
    #Discretized solution parameters
    numerical.N      = 101
    numerical.pie    = N_Pe/numerical.N #Block Peclet number (u*dx)/(D*phi)
    numerical.cr     = 0.5              #Dimensionless Courant number
    numerical.eta    = numerical.cr/numerical.pie #Dimensionless block size
    numerical.dxD    = 1/numerical.N    #Dimensionless block size
    numerical.PVI    = 2.0
    numerical.xD     = np.empty((numerical.N,1))

    numerical.dx     = 1/numerical.N
    numerical.dt     = 2.1064e3*3/(numerical.N**2 / 25)
    
    #position of the block centres
    for i in range(0,numerical.N):
        numerical.xD[i,0] = ((1.0/numerical.N)*(i+1.0)) - (1.0/numerical.N)/2
        
    numerical.method = 'explicit'  #'Implicit' or 'explicit'
    numerical.scheme = 'upwind'    #'upwind' or 'center'
    
    #Initial condition
    IC.C = np.zeros((numerical.N, 1))  #Concentration

    #Wells
    # Type of wells
    well_type = 4  #1= left to right, 2=right to left, 3=injector in center, 4=producer in center
    IC.P = np.empty((numerical.N, 1))
    TAO = lil_matrix((numerical.N, numerical.N))
    
    if well_type==1:
        TAO[0,0] = 1.0
        TAO[numerical.N-1,numerical.N-1] = -1.0  
        for i in range(0,numerical.N): 
            IC.P[i,0] = numerical.N - i
            
    elif well_type==2:
        TAO[0,0] = -1.0
        TAO[numerical.N-1,numerical.N-1] = 1.0  
        for i in range(0,numerical.N): 
            IC.P[i,0] = i - numerical.N  
            
    elif well_type==3:
        TAO[0,0] = -1.0
        TAO[numerical.N-1,numerical.N-1] = -1.0 
        TAO[(numerical.N-1)/2,(numerical.N-1)/2] =  2.0  
        for i in range(0,numerical.N-1): 
            IC.P[i,0] = -(i-(numerical.N-1)/2)**2.0 
        IC.P[numerical.N-1,0] =  IC.P[0,0]
    
    elif well_type==4:
        TAO[0,0] = 1.0
        TAO[numerical.N-1,numerical.N-1] = 1.0 
        TAO[(numerical.N-1)/2,(numerical.N-1)/2] =-2.0  
        for i in range(0,numerical.N-1): 
            IC.P[i,0] = (i-(numerical.N-1)/2)**2.0 
        IC.P[numerical.N-1,0] =  IC.P[0,0]

    return numerical, reservoir, fluid, IC, N_Pe, TAO, well_type