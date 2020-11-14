"""
reservoir simulation assignment 2
1D reservoir simulation Q3.4: Main file
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""

#import inbuilt libraries
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import time as timer

#importing personal libraries
from input_file import inputfile
from post_process import post_process #, analytical

#making fluid, reservoir,simulation, grid, BC and IC classes
class fluid:
    def __init__(self):
        self.mu = []
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

inputfile(fluid,reservoir,numerical,BC,IC)     #uploading all input properties          

#Setting up matrix T, B, and Q
T = lil_matrix((numerical.N, numerical.N))
B = lil_matrix((numerical.N, numerical.N))
Q = lil_matrix((numerical.N, 1))

trans = reservoir.k * reservoir.A /(fluid.mu * fluid.Bw * numerical.dx) #Transmissibility vector 

start = timer.clock()

for i in range(0,numerical.N):
    if i==0:               #left boundary
        T[i,i+1] = - trans[i]
        if 'Neumann' in BC.type[0]:
            T[i,i] = T[i,i] - T[i,i+1]
        elif 'Dirichlet' in BC.type[0]:
            T[i,i] = T[i,i] - T[i,i+1] + 2*trans[i]
            Q[i,0]   = 2 * trans[i] * (BC.value[0][0] - (fluid.rho/144.0)*BC.depth[0][0] )* 6.33E-3      
    
    elif i==numerical.N-1:  #right boundary
        T[i,i-1] = - trans[i]
        if 'Neumann' in BC.type[1]:
            T[i,i] = T[i,i] - T[i,i-1]
        elif 'Dirichlet' in BC.type[1]:
            T[i,i] = T[i,i] - T[i,i-1] + 2*trans[i]
            Q[i,0]   = 2 * trans[i] * (BC.value[1][0] - (fluid.rho/144.0)*BC.depth[1][0] ) * 6.33E-3      

    else:
         T[i,i-1] = -trans[i]
         T[i,i+1] = -trans[i]                       
         T[i,i]   = -(T[i,i-1] + T[i,i+1])
     
    B[i,i] = reservoir.A * numerical.dx[i,0] * reservoir.phi[i,0] * fluid.ct / fluid.Bw #accumulation

T = (6.33E-3 * T).tocsr() #multiplying with the conversion factor 
G = csr_matrix((fluid.rho/144.0)*(T@numerical.D))

# time marching initialization
time = 0     #initializing time
n    = 0     #time step initialization
nmax = int(numerical.tfinal / numerical.dt)   #total number of time steps
time = np.zeros((nmax+1))                     #initializing time vector
P_plot      = np.zeros((numerical.N, nmax+1)) #matrix to save pressure 
P           = (IC.P).copy()                   #initializing iterable current pressure
P_plot[:,n] = P[:,0]                          #saving initial pressure
B    = B.tocsr()
Q    = Q.tocsr()

# time marching
for i in range(n,nmax):
    P_old = csr_matrix(P) #dummy vector as a placeholder for old P  
    IM = (1-numerical.theta)*T + B / numerical.dt             #implicit part coefficient in Eq. 3.44   
    EX = ((B / numerical.dt - numerical.theta*T) @ P_old)+Q+G #explicit part or RHS of Eq. 3.44
    
    P = np.transpose([spsolve(IM,EX)])                        #solving IM*P = EX or Ax=B     
    
    P_plot[:,i+1] = P[:,0]                                    #storing the pressure values inside the domain
    time[i+1] = time[i] + numerical.dt                        #marching time

end = timer.clock()

print(end-start)

post_process(P_plot,numerical.xc[:,0],time,10)