"""
reservoir simulation assignment 3
1D reservoir simulation: Arrays
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/29/2020
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
def myarrays3(reservoir,fluid,numerical,IC):

    #Setting up matrix T, B, and Q
    ETA = lil_matrix((numerical.N, numerical.N))
    CR  = lil_matrix((numerical.N, numerical.N))
    Q   = lil_matrix((numerical.N, 1))
    C   = lil_matrix((numerical.N, 1))

    for i in range(0,numerical.N):
        if i == 0:                 #left boundary  
            ETA[i,i+1] = -1
            ETA[i,i]   =  1    

            if IC.P[i,0] > IC.P[i+1,0]:
                CR[i,i]  = CR[i,i] + 1.0*courant(i,i+1,IC.P,numerical,fluid,reservoir)
            elif IC.P[i,0] <= IC.P[i+1,0]:
                CR[i,i+1]= CR[i,i+1] - 1.0*courant(i,i+1,IC.P,numerical,fluid,reservoir)

        elif i==numerical.N-1:     #right boundary
            ETA[i,i-1] = -1
            ETA[i,i]   =  1    
            
            if IC.P[i,0] < IC.P[i-1,0]:
                CR[i,i-1]  = -1.0*courant(i,i-1,IC.P,numerical,fluid,reservoir)
            elif IC.P[i,0] >= IC.P[i-1,0]:
                CR[i,i]    =   CR[i,i] + 1.0*courant(i,i-1,IC.P,numerical,fluid,reservoir)
       

        else:                       #interior
            
            if IC.P[i,0] < IC.P[i-1,0]:
                CR[i,i-1]  = -1.0*courant(i,i-1,IC.P,numerical,fluid,reservoir)
            else:
                CR[i,i]    =  CR[i,i] + 1.0*courant(i,i-1,IC.P,numerical,fluid,reservoir)
        
            if IC.P[i,0] > IC.P[i+1,0]:
                CR[i,i]  = CR[i,i] + 1.0*courant(i,i+1,IC.P,numerical,fluid,reservoir)
            else:
                CR[i,i+1]  = CR[i,i+1] - 1.0*courant(i,i+1,IC.P,numerical,fluid,reservoir)   
            
            ETA[i,i+1] = -1
            ETA[i,i-1] = -1
            ETA[i,i]   =  2   
        
    return ETA,CR;  

numerical.N = 5               

# Type of wells
well_type = 4  #1= left to right, 2=right to left, 3=injector in center, 4=producer in center
IC.P = np.empty((numerical.N, 1))
TAO = lil_matrix((numerical.N, numerical.N))

fluid.visc = 1.0
reservoir.perm = 1.0
numerical.dx = 1/numerical.N

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

ETA,CR = myarrays3(reservoir,fluid,numerical,IC) 