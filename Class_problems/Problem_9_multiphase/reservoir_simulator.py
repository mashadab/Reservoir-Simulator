"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Main file
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

#import inbuilt libraries
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import time as timer

#importing personal libraries
from input_file import inputfile
from myarrays import myarrays  #importing the arrays
from init_plot import initial_plot #for plotting
from updatewells import updatewells #for updating the wells
from postprocess import postprocess #for plotting wells

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
class well:
    def __init__(self):
        self.xmin = []
        self.block = []

start = timer.clock()

inputfile(fluid,reservoir,numerical,BC,IC,well)  #uploading all input properties          

#Calculating the arrays
[T, B, Q, G] = myarrays(fluid,reservoir,numerical,BC,IC)

#Updating the well
[J, Q]       = updatewells(reservoir,fluid,numerical,IC.P,well,Q)

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
well.typetime       = np.kron(np.ones((nmax,1)),np.transpose(well.type))
well.constrainttime = np.kron(np.ones((nmax,1)),np.transpose(well.constraint))

# time marching
for i in range(n,nmax):
    P_old = csr_matrix(P) #dummy vector as a placeholder for old P

    #Switching according to the problem
    ##############################################################################
    #Switching after 500 days
    if time[i] > numerical.tswitch:
        well.type= [[2],[2],[2],[2],[2],[2],[2],[2],[1],[1]] # 1 for rate, 2 for BHP
        well.constraint = [[502.5],[502.5],[502.5], [502.5], [502.5], [502.5], [502.5], [502.5],[500*5.61], [500*5.61]]  # rate = scf/day (+ for injector); BHP = psi (always +ve)
    else:
        # Well 1 BHP
        if well.constraint[0][0] / well.Jvec[0][0] + P[well.block[0][0],0] <= 502.5:
        #if P[well.block[0][0],0] <= 502.5:
            well.type[0][0] = 2
            well.constraint[0][0] = 502.5

        #Well 2 BHP
        if well.constraint[1][0] / well.Jvec[1][0] + P[well.block[1][0],0] <= 502.5:
        #if P[well.block[1][0],0] <= 502.5:
            well.type[1][0] = 2
            well.constraint[1][0] = 502.5

    well.typetime[i,:] = np.copy(np.transpose(well.type[:][:]))
    well.constrainttime[i,:] = np.copy(np.transpose(well.constraint[:][:]))

    #Updating the well
    [J, Q]       = updatewells(reservoir,fluid,numerical,IC.P,well,BC)
    ##############################################################################

    IM = (1-numerical.theta)*(T + J) + B / numerical.dt             #implicit part coefficient in Eq. 3.44   
    EX = ((B / numerical.dt - numerical.theta*(T + J)) @ P_old)+Q+G # explicit part or RHS of Eq. 3.44
    
    P = np.transpose([spsolve(IM,EX)])                        #solving IM*P = EX or Ax=B     
    
    P_plot[:,i+1] = P[:,0]                                    #storing the pressure values inside the domain
    time[i+1] = time[i] + numerical.dt                        #marching time

end = timer.clock()
print(end-start)

#Create the plots
initial_plot(reservoir,numerical,P_plot[:,0],time[0])
initial_plot(reservoir,numerical,P_plot[:,300],time[300])
initial_plot(reservoir,numerical,P_plot[:,600],time[600])
initial_plot(reservoir,numerical,P_plot[:,nmax],time[nmax])

postprocess(P_plot,numerical,well,time)