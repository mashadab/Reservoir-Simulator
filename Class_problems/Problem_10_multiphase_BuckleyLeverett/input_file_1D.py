"""
reservoir simulation project 2 (2020): Problem 1
2D Multiphase reservoir simulation: Input file
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 11/12/2020
"""

import numpy as np
import matplotlib.pyplot as plt

from rel_perm import rel_perm
from Thalf import Thalf
from myarrays import myarrays

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

class petro:
    def __init__(self):
        self.xmin = []        

class well:
    def __init__(self):
        self.xmin = []
        self.xblock = []

#fluid, reservoir and simulation parameters   
def inputfile(fluid,reservoir,petro,numerical,BC,IC,well):
    # Numerical simulation parameters
    numerical.dt     = 10   #time step (days)
    numerical.tfinal = 4    #final time [days]
    numerical.PV_final=1.0  #non-dimensional pore volume injected (a non dimensional time)
    numerical.Nx     = 50   #number of grid blocks in x-direction
    numerical.Ny     = 1    #number of grid blocks in y-direction 
    numerical.N  = numerical.Nx * numerical.Ny #Total number of grid blocks
    #numerical.theta  = 1.0   #type of method theta = 1-Explicit, 0-Implicit, 0.5-CN
    numerical.method = 'IMPES' #Implicit pressure explicit saturation solver

    # Fluid parameters
    fluid.muw = 1.0*np.ones((numerical.N, 1))   #fluid viscosity [centipoise]  
    fluid.Bw  = 1.0*np.ones((numerical.N, 1))   #formation volume factor of water [rb/stb]  
    fluid.cw  = 1.0E-5                          #total compressibility: rock + fluid [1/psi] 

    fluid.muo = 1.0*np.ones((numerical.N, 1)) #fluid viscosity [centipoise]  
    fluid.Bo  = 1.0*np.ones((numerical.N, 1)) #formation volume factor of oil [rb/stb]  
    fluid.co  = 1.0E-5     #total compressibility: rock + fluid [1/psi] 

    fluid.rho = 62.4        #density of the fluid [lbm/ft^3]
    fluid.relperm = 1.0*np.ones((numerical.N, 1))  #Relative permeability of the fluid

    # Multiphase/Relative permeability values and Capillary pressure
    petro.Swr  = 0.2    #residual water saturation
    petro.Sor  = 0.2    #residual oil saturation
    petro.nw   = 3.0    #Corey-Brooks exponent (water)
    petro.no   = 3.0    #Corey-Brooks exponent (oil)
    petro.krw0 = 0.2    #Corey-Brooks endpoint (water)
    petro.kro0 = 1.0    #Corey-Brooks endpoint (oil)
    petro.lamda= 2.0    #fitting parameter for Corey-Brooks model
    petro.Pe   = 3.5    #capillary entry pressure [psi]
    
    '''
    #reading files
    depth    =-np.loadtxt("PJ1-Depth.txt")
    porosity = np.loadtxt("PJ1-Porosity.txt")
    permx    = np.loadtxt("PJ1-Permeability.txt")    
    '''                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    
    # Reservoir parameters    
    reservoir.L  = 1000.0  #length of the reservoir [feet] 
    reservoir.h  = 10.0    #height of the reservoir [feet] 
    reservoir.W  = 1000.0  #width of the reservoir [feet]    
    reservoir.phi   = 0.2*np.ones((numerical.N,1))  #porosity of the reservior vector [mDarcy] 
    reservoir.permx = 100.0*np.ones((numerical.N,1))#fluid permeability in x direction vector [mDarcy] 
    reservoir.permx[reservoir.permx <= 1E-6] = 1E-16
    reservoir.permy  = 1.0*reservoir.permx #fluid permeability in y direction vector [mDarcy]                      
    reservoir.permz  = 1.0*reservoir.permx #fluid permeability in z direction vector [mDarcy]   
    reservoir.Dref   = 7495              #reference depth [feet]
    reservoir.alpha  = 0.0*np.pi/6.0     #dip angle [in radians]
    reservoir.cfr    = 0.0               #formation of rock compressibility    
                                         
    # Well parameters
    well.x   = [[0.001*reservoir.L],[0.999*reservoir.L]]
    well.y   = [[reservoir.W/2],[reservoir.W/2]]    
    well.type= [[1],[1]] # 1 for rate, 2 for BHP
    well.constraint = [[426.5],[-426.5]]   # rate = scf/day (+ for injector); BHP = psi (always +ve)
    well.rw  = [[0.25],[0.25]]             # well radius, ft
    well.skin= [[0],[0]] # well skin friction factor, dimensionless
    well.direction = [['v'],['v']]   # direction of the well: vertical v or horizontal hx or hy

    #Defining numerical parameters for discretized solution
    numerical.dx1 = (reservoir.L/numerical.Nx)*np.ones((numerical.Nx, 1)) #block thickness in x vector
    numerical.dy1 = (reservoir.W/numerical.Ny)*np.ones((numerical.Ny, 1)) #block thickness in y vector
    [numerical.dX,numerical.dY] = np.meshgrid(numerical.dx1,numerical.dy1) 
    
    numerical.dx  = np.reshape(numerical.dX, (numerical.N,1))    #building the single dx column vector
    numerical.dy  = np.reshape(numerical.dY, (numerical.N,1))    #building the single dy column vector

    #position of the block centres x-direction
    numerical.xc = np.empty((numerical.Nx, 1))
    numerical.xc[0,0] = 0.5 * numerical.dx[0,0]
    for i in range(1,numerical.Nx):
        numerical.xc[i,0] = numerical.xc[i-1,0] + 0.5*(numerical.dx1[i-1,0] + numerical.dx1[i,0])

    #position of the block centres y-direction
    numerical.yc = np.empty((numerical.Ny, 1))
    numerical.yc[0,0] = 0.5 * numerical.dy[0,0]
    for i in range(1,numerical.Ny):
        numerical.yc[i,0] = numerical.yc[i-1,0] + 0.5*(numerical.dy1[i-1,0] + numerical.dy1[i,0])
    
    [numerical.Xc,numerical.Yc] = np.meshgrid(numerical.xc,numerical.yc) 
    
    numerical.x1  = np.reshape(numerical.Xc, (numerical.N,1))    #building the single X column vector
    numerical.y1  = np.reshape(numerical.Yc, (numerical.N,1))    #building the single Y column vector

    #depth vector
    numerical.D = np.array([[4300],[4300],[4300]])

    BC.type  = [['Neumann'],['Neumann'],['Neumann'],['Neumann']] #type of BC: left, right, bottom, top
    BC.value = [[0],[0],[0],[0]]    #value of the boundary condition: psi or ft^3/day
    #BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)],[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    IC.P     = 1000.0*np.ones((numerical.N,1)) #Initial Pressure
    IC.Sw    = 0.2*np.ones((numerical.N,1))  #Initial saturation
                                       
#depth    =-np.loadtxt("PJ1-Depth.txt")
#porosity = np.loadtxt("PJ1-Porosity.txt")
#permx    = np.loadtxt("PJ1-Permeability.txt") 

inputfile(fluid,reservoir,petro,numerical,BC,IC,well)

print(rel_perm(petro,0.5))

print(Thalf(0,1,'x',fluid,reservoir,petro,numerical,IC.P,IC.Sw))

Tw, To, T, d11, d12, d21, d22, D, G = myarrays(fluid,reservoir,petro,numerical,BC,IC.P,IC.Sw)