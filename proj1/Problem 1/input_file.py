"""
reservoir simulation project 1 (2020): Problem 1
2D reservoir simulation: Input file
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
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
    reservoir.L  = 1200.0 #length of the reservoir [feet] 
    reservoir.h  = 200.0  #height of the reservoir [feet] 
    reservoir.W  = 600.0  #width of the reservoir [feet]    
    reservoir.phi= np.transpose(np.array([[0.26, 0.20, 0.23, 0.22, 0.24, 0.18, 0.25, 0.20, 0.22]]))     #porosity of the reservior vector
    reservoir.permx  = np.transpose(np.array([[1800, 2000, 1600, 2000, 2500, 1200, 1000, 2000, 2200]])) #fluid permeability in x direction vector [mDarcy] 
    reservoir.permy  = 2.0*reservoir.permx     
    #reservoir.permx  = np.array([[50],[40],[30],[60],[50],[40],[70],[60],[50]],dtype='float64')  #fluid permeability in x direction vector [mDarcy] 
    #reservoir.permy  = 2.0*reservoir.permx                                                       #fluid permeability in y direction vector [mDarcy]                      
    reservoir.Dref   = 1000/(fluid.rho/144.0)     #reference depth [feet]
    reservoir.alpha  = 0.0*np.pi/6.0              #dip angle [in radians]

    # Well parameters
    well.x  = [[600.0],[1000.0]]
    well.y  = [[300.0],[500.0]]    
    well.type= [[1],[2]] # 1 for rate, 2 for BHP
    well.constraint = [[-10000],[1500]]  # rate = scf/day (+ for injector); BHP = psi (always +ve)
    well.rw = [[0.5],[0.5]]             # well radius, ft
    well.skin =[[-0.75],[-0.75]]        # well skin friction factor, dimensionless

    #Defining numerical parameters for discretized solution
    numerical.dx1  = np.transpose(np.array([[200, 600, 400]]))    #building the dx  vector
    numerical.dy1  = np.transpose(np.array([[100, 300, 200]]))    #building the dy  vector

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
    numerical.D = np.transpose(np.array([[1883.2, 2275.6, 2300.6, 1883.2, 2275.6, 2300.6, 1883.2, 2275.6, 2300.6]]))
    
    BC.type  = [['Neumann'],['Dirichlet'],['Neumann'],['Neumann']] #type of BC: left, right, bottom, top
    BC.value = [[0],[1200],[0],[0]]    #value of the boundary condition: psi or ft^3/day
    BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)],[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    #IC.P = 1000.0*np.ones((numerical.N, 1)) #initial condition: pressure in the domain
    IC.P = 1000.0 + (fluid.rho/144.0)*(numerical.D - reservoir.Dref)