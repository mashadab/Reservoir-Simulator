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
class well:
    def __init__(self):
        self.xmin = []
        self.xblock = []

#fluid, reservoir and simulation parameters   
def inputfile(fluid,reservoir,numerical,BC,IC,well):
    # Numerical simulation parameters
    numerical.dt     = 0.01   #time step (days)
    numerical.tfinal = 20     #final time [days]
    numerical.Nx     = 3     #number of grid blocks in x-direction
    numerical.Ny     = 3     #number of grid blocks in y-direction    
    numerical.N  = numerical.Nx * numerical.Ny #Total number of grid blocks
    numerical.theta  =  0.0      #type of method theta = 1-Explicit, 0-Implicit, 0.5-CN

    # Fluid parameters
    fluid.mu  = 1.0*np.ones((numerical.N, 1)) #fluid viscosity [centipoise]  
    fluid.Bw  = 1.0*np.ones((numerical.N, 1)) #formation volume factor of water [rb/stb]  
    fluid.ct  = 5E-6        #total compressibility [1/psi] 
    fluid.rho = 62.4        #density of the fluid [lbm/ft^3]
    fluid.relperm = 1.0*np.ones((numerical.N, 1))  #Relative permeability of the fluid

    # Reservoir parameters    
    reservoir.L  = 300.0  #length of the reservoir [feet] 
    reservoir.h  = 100.0  #height of the reservoir [feet] 
    reservoir.W  = 300.0  #width of the reservoir [feet]    
    reservoir.phi= 0.2*np.ones((numerical.N, 1))     #porosity of the reservior vector
    reservoir.permx  = 100*np.ones((numerical.N, 1)) #fluid permeability in x direction vector [mDarcy] 
    reservoir.permy  = np.copy(reservoir.permx)      
    #reservoir.permx  = np.array([[50],[40],[30],[60],[50],[40],[70],[60],[50]],dtype='float64')  #fluid permeability in x direction vector [mDarcy] 
    #reservoir.permy  = 2.0*reservoir.permx                                                       #fluid permeability in y direction vector [mDarcy]                      
    reservoir.Dref   = 1000/(fluid.rho/144.0)     #reference depth [feet]
    reservoir.alpha  = 0.0*np.pi/6.0              #dip angle [in radians]

    # Well parameters
    well.x  = [[reservoir.L/2.0],[299.0]]
    well.y  = [[reservoir.W/2.0],[299.0]]    
    well.type= [[1],[2]] # 1 for rate, 2 for BHP
    well.constraint = [[-1000],[1500]]  # rate = scf/day (+ for injector); BHP = psi (always +ve)
    well.rw = [[0.5],[0.5]]             # well radius, ft
    well.skin =[[-0.75],[-0.75]]        # well skin friction factor, dimensionless

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
        numerical.xc[i,0] = numerical.xc[i-1,0] + 0.5*(numerical.dx[i-1,0] + numerical.dx[i,0])

    #position of the block centres y-direction
    numerical.yc = np.empty((numerical.Ny, 1))
    numerical.yc[0,0] = 0.5 * numerical.dy[0,0]
    for i in range(1,numerical.Ny):
        numerical.yc[i,0] = numerical.xc[i-1,0] + 0.5*(numerical.dy[i-1,0] + numerical.dy[i,0])
    
    [numerical.Xc,numerical.Yc] = np.meshgrid(numerical.xc,numerical.yc) 
    
    numerical.x1  = np.reshape(numerical.Xc, (numerical.N,1))    #building the single X column vector
    numerical.y1  = np.reshape(numerical.Yc, (numerical.N,1))    #building the single Y column vector
    
    #depth vector 
    numerical.D = reservoir.Dref-numerical.x1*np.sin(reservoir.alpha)
    
    BC.type  = [['Neumann'],['Dirichlet'],['Neumann'],['Neumann']] #type of BC: left, right, bottom, top
    BC.value = [[0],[1200],[0],[0]]    #value of the boundary condition: psi or ft^3/day
    BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)],[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    IC.P = 1000.0*np.ones((numerical.N, 1)) #initial condition: pressure in the domain
    #IC.P = 1000.0 + (fluid.rho/144.0)*(numerical.D - reservoir.Dref)