"""
reservoir simulation assignment 2
1D reservoir simulation Q3.4: Input files
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""
import numpy as np

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

#fluid, reservoir and simulation parameters   
def inputfile(fluid,reservoir,numerical,BC,IC):
    fluid.mu  = 1.0         #fluid viscosity [centipoise]  
    fluid.Bw  = 1.0         #formation volume factor of water [rb/stb]  
    fluid.ct  = 1E-6        #total compressibility [1/psi] 
    fluid.rho = 62.4        #density of the fluid [lbm/ft^3]
    
    numerical.dt     = 1.0    #time step (days)
    numerical.tfinal = 200     #final time [days]
    numerical.N      = 40000    #number of grid blocks
    numerical.theta  =  0.0      #type of method theta = 1-Explicit, 0-Implicit, 0.5-CN
    
    reservoir.L  = 10000.0  #length of the reservoir [feet] 
    reservoir.A  = 2E5      #cross sectinal area of the reservoir [feet^2]     
    reservoir.phi= 0.2*np.ones((numerical.N, 1))   #porosity of the reservior vector
    reservoir.k  = 50.0*np.ones((numerical.N, 1))  #fluid permeability vector [mDarcy] 
    reservoir.Dref=1000/(fluid.rho/144.0)      #reference depth [feet]
    reservoir.alpha=1.0*np.pi/6.0              #dip angle [in radians]
    
    numerical.dx = (reservoir.L/numerical.N)*np.ones((numerical.N, 1)) #block thickness vector
    numerical.xc = np.empty((numerical.N, 1))
    numerical.xc[0,0] = 0.5 * numerical.dx[0,0]
    
    #position of the block centres
    for i in range(0,numerical.N):
        numerical.xc[i,0] = numerical.xc[i-1,0] + 0.5*(numerical.dx[i-1,0] + numerical.dx[i,0])
    
    #depth vector 
    numerical.D = reservoir.Dref-numerical.xc*np.sin(reservoir.alpha)
    
    BC.type  = [['Dirichlet'],['Neumann']]  #type of BC as first cell, last cell
    BC.value = [[2000],[0]]                 #value of the boundary condition: psi or ft^3/day
    BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    IC.P = 1000.0*np.ones((numerical.N, 1)) #initial condition: pressure in the domain
    #IC.P = 1000.0 + (fluid.rho/144.0)*(numerical.D - reservoir.Dref)