"""
reservoir simulation assignment 3
1D reservoir simulation Q7: Main file (Implicit pressure explicit concentration)
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/29/2020
"""

#import inbuilt libraries
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, identity
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import time as timer
from math import floor, ceil

#importing personal libraries
from input_file1D_new import inputfile
from myarray4 import myarrays4

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

numerical, reservoir, fluid, IC, N_Pe, TAO, well, well_type = inputfile()     #uploading all input properties          

#creating arrays
ETA,CR = myarrays4(reservoir,fluid,numerical,IC)
I      = identity(numerical.N)

#looping through time
C_star = np.zeros((numerical.N,1))

Nexit  = floor(1 + (numerical.PVI*numerical.N/numerical.cr))
t_D = np.empty((100000))
t_D[0]= 0
k   = 0
Cexit = np.zeros((100000,1))
k1  = 1.0
C   = np.copy(IC.C)

C_plot= np.zeros((numerical.N, 100000)) #matrix to save pressure 
C_plot[:,0] = IC.C[:,0] 

while (t_D[k] < numerical.PVI):      
    if well_type==1:     #left to right
        C_star[0,0] = 1.0
        C_star[numerical.N-1,0] = C[numerical.N-1,0]   

    elif well_type==2:     #right to left
        C_star[0,0] = C[0,0] 
        C_star[numerical.N-1,0] = 1.0  

    elif well_type==3:     #injector in the middle
        C_star[0,0] = C[0,0] 
        C_star[numerical.N-1,0] = C[numerical.N-1,0] 
        C_star[ceil((numerical.N-1)/2),0] = 1.0

    elif well_type==4:     #producer in the middle
        C_star[0,0] = 1.0 
        C_star[numerical.N-1,0] = 1.0
        C_star[ceil((numerical.N-1)/2),0] = C[ceil((numerical.N-1)/2),0]

    if numerical.method == 'explicit':
        C =  C - (CR + ETA) @ C + TAO @ C_star

    k = k+1
    C_plot[:,k] = C[:,0]   
    t_D[k]= CR[0,0]*(k-1)/numerical.N
    Cexit[k] = C[numerical.N-1,0]

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
plt.plot(t_D[0:k-1],Cexit[0:k-1])
plt.xlabel(r'$t_D$')
plt.ylabel(r'$C_D$')
plt.savefig('CdvsTd.png',bbox_inches='tight', dpi = 600)