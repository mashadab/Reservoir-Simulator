"""
reservoir simulation assignment 2
1D reservoir simulation Q3.4: Main file
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""

#import inbuilt libraries
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt

#importing personal libraries
from input_file import inputfile
from post_process import post_process

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
T = sp.csr_matrix((numerical.N, numerical.N), dtype=np.float64)
B = sp.csr_matrix((numerical.N, numerical.N), dtype=np.float64)
Q = sp.csr_matrix((numerical.N, 1), dtype=np.float64)

trans = reservoir.k * reservoir.A /(fluid.mu * fluid.Bw * numerical.dx) #Transmissibility vector 

for i in range(0,numerical.N):
    if i==0:               #left boundary
        T[i,i+1] = - trans[i]
        if 'Neumann' in BC.type[0]:
            T[i,i] = T[i,i] - T[i,i+1]
        elif 'Dirichlet' in BC.type[0]:
            T[i,i] = T[i,i] - T[i,i+1] + 2*trans[i]
            Q[i,0]   = 2 * trans[i] * (BC.value[0]) * 6.33E-3      
    
    elif i==numerical.N-1:  #right boundary
        T[i,i-1] = - trans[i]
        if 'Neumann' in BC.type[1]:
            T[i,i] = T[i,i] - T[i,i-1]
        elif 'Dirichlet' in BC.type[1]:
            T[i,i] = T[i,i] - T[i,i-1] + 2*trans[i]
            Q[i,0]   = 2 * trans[i] * BC.value[1] * 6.33E-3      

    else:
         T[i,i-1] = -trans[i]
         T[i,i+1] = -trans[i]                       
         T[i,i]   = -(T[i,i-1] + T[i,i+1])
     
    B[i,i] = reservoir.A * numerical.dx[i,0] * reservoir.phi[i,0] * fluid.ct / fluid.Bw #accumulation

T = 6.33E-3 * T #multiplying with the conversion factor 

# time marching initialization
time = 0     #initializing time
n    = 0     #time step initialization
nmax = int(numerical.tfinal / numerical.dt) #total number of time steps
time = np.zeros((nmax+1))                     #initializing time vector
P_plot      = np.zeros((numerical.N, nmax+1)) #matrix to save pressure 
P           = (IC.P).copy()                   #initializing iterable current pressure
P_plot[:,n] = P[:,0]                          #saving initial pressure

# time marching
for i in range(n,nmax):
    P_old       = sp.csr_matrix(P) #dummy vector as a placeholder for old P
    
    IM = (1-numerical.theta)*T + B / numerical.dt             #implicit part coefficient in Eq. 3.44
    EX = ((B / numerical.dt - numerical.theta*T) @ P_old) + Q #explicit part or RHS of Eq. 3.44
    
    P = np.transpose([linalg.spsolve(IM,EX)])                 #solving IM*P = EX or Ax=B
    P_plot[:,i+1] = P[:,0]                                    #storing the pressure values inside the domain
    time[i+1] = time[i] + numerical.dt                        #marching time

post_process(P_plot,numerical.xc[:,0],time,nmax)

'''
#plotting
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(numerical.xc[:,0],P_plot[:,0],'r-',label=r"$%0.1f$ [days]" % time[0])
plt.plot(numerical.xc[:,0],P_plot[:,10],'r--',label=r"$%0.1f$ [days]" % time[10])
plt.plot(numerical.xc[:,0],P_plot[:,20],'r-.',label=r"$%0.1f$ [days]" % time[20])
plt.plot(numerical.xc[:,0],P_plot[:,30],'b-',label=r"$%0.1f$ [days]" % time[30])
plt.plot(numerical.xc[:,0],P_plot[:,40],'b--',label=r"$%0.1f$ [days]" % time[40])
plt.plot(numerical.xc[:,0],P_plot[:,50],'b-.',label=r"$%0.1f$ [days]" % time[50])
plt.plot(numerical.xc[:,0],P_plot[:,100],'g-',label=r"$%0.1f$ [days]" % time[100])
plt.plot(numerical.xc[:,0],P_plot[:,150],'g--',label=r"$%0.1f$ [days]" % time[150])
plt.plot(numerical.xc[:,0],P_plot[:,200],'g-.',label=r"$%0.1f$ [days]" % time[200])
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$Pressure$ [psi]')
plt.xlabel(r'$x$ [feet]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig('Pvsx.png',bbox_inches='tight', dpi = 600)
'''