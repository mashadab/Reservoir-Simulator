"""
reservoir simulation assignment 3
2D reservoir simulation Q7: Main file (Implicit pressure explicit saturation)
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 11/17/2020
"""

#import inbuilt libraries
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, identity
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import time as timer
from math import floor, ceil

#importing personal libraries
from input_file import inputfile
from myarrays import myarrays
from updatewells import updatewells


#making simulation and IC classes
class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = [] 
class fluid:
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
class petro:
    def __init__(self):
        self.xmin = []
class well:
    def __init__(self):
        self.xmin = []

#loading inputfile
inputfile(fluid,reservoir,petro,numerical,BC,IC,well)
#Implicit pressure and explicit saturation for update
#looping through time
t = np.empty((100000))
t[0]= 0
k   = 0
P   = np.copy(IC.P)
Sw  = np.array(np.copy(IC.Sw))
nmax= ceil(numerical.tfinal / numerical.dt)

P_plot= np.zeros((numerical.N, nmax + 1)) #matrix to save pressure 
P_plot[:,0] = IC.P[:,0] 
Sw_plot= np.zeros((numerical.N, nmax + 1)) #matrix to save pressure 
Sw_plot[:,0]= IC.Sw[:,0] 

while (t[k] < numerical.tfinal):      
    P_old = np.copy(P)   #Placeholdering the old array
    
    #Calculating the arrays
    Tw, To, T, d11, d12, d21, d22, D, G = myarrays(fluid,reservoir,petro,numerical,BC,P,Sw)
    
    #updating the wells
    well, Qw, Qo, J = updatewells(reservoir,fluid,numerical,petro,P,Sw,well)
    Q = (-d22 @ inv(d12)) @ Qw + Qo
    
    if numerical.method == 'IMPES':
        IM = T + D              #implicit part coefficient in Eq. 3.44   
        EX = D @ P_old + Q + G  #explicit part or RHS of Eq. 3.44
    
        P = np.transpose([spsolve(IM,EX)])                              #solving IM*P = EX or Ax=B        
        Sw = Sw + inv(d12) @ (-Tw @ P - d11 @ (P - P_old) + Qw)         #explicit saturation

    k = k+1
    P_plot[:,k] = P[:,0]  
    Sw_plot[:,k]= np.array(Sw)[:,0]
    t[k]= t[k-1] + numerical.dt


P_plot[np.argwhere(reservoir.permx < 0.01)] = np.nan

#post process
'''
plt.figure()
plt.plot(numerical.xD,C_plot[:,21],'r-',label=f'Time=%0.2f' % t_D[21])
plt.plot(numerical.xD,C_plot[:,102],'b-',label=f'Time=%0.2f' % t_D[102])
plt.plot(numerical.xD,C_plot[:,203],'k-',label=f'Time=%0.2f' % t_D[203])
plt.legend(loc='best', shadow=False, fontsize='medium')
plt.xlabel(r'$x_D$')
plt.ylabel(r'$C_D$')
'''

plt.figure()
plt.plot(numerical.xc,P_plot)
plt.xlabel(r'$x$ [feet]')
plt.ylabel(r'$P$ [psi]')
plt.savefig('PvsT.png',bbox_inches='tight', dpi = 600)