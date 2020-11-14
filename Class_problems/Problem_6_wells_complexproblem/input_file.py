"""
reservoir simulation assignment 3
1D reservoir simulation Q3.4: Input files
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""
import numpy as np
import matplotlib.pyplot as plt

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
    numerical.dt     = 1     #time step (days)
    numerical.tfinal = 1000  #final time [days]
    numerical.Nx     = 80    #number of grid blocks in x-direction
    numerical.Ny     = 75    #number of grid blocks in y-direction 
    numerical.N  = numerical.Nx * numerical.Ny #Total number of grid blocks
    numerical.theta  = 0.0   #type of method theta = 1-Explicit, 0-Implicit, 0.5-CN
    numerical.tswitch= 500   #time after which the wells switch [days]  

    # Fluid parameters
    fluid.mu  = 0.383*np.ones((numerical.N, 1)) #fluid viscosity [centipoise]  
    fluid.Bw  = 1.012298811*np.ones((numerical.N, 1))#formation volume factor of water [rb/stb]  
    fluid.ct  = 3.87E-6     #total compressibility: rock + fluid [1/psi] 
    fluid.rho = 62.4        #density of the fluid [lbm/ft^3]
    fluid.relperm = 1.0*np.ones((numerical.N, 1))  #Relative permeability of the fluid
    
    #reading files
    depth    =-np.loadtxt("PJ1-Depth.txt")
    porosity = np.loadtxt("PJ1-Porosity.txt")
    permx    = np.loadtxt("PJ1-Permeability.txt")    

    # Reservoir parameters    
    reservoir.L  = 6000.0  #length of the reservoir [feet] 
    reservoir.h  = 50.0    #height of the reservoir [feet] 
    reservoir.W  = 7500.0  #width of the reservoir [feet]    
    reservoir.phi   = np.reshape(porosity, (numerical.N,1))  #porosity of the reservior vector [mDarcy] 
    reservoir.permx = np.reshape(permx, (numerical.N,1))  #fluid permeability in x direction vector [mDarcy] 
    reservoir.permx[reservoir.permx <= 1E-6] = 1E-16
    reservoir.permy  = 0.15*reservoir.permx               #fluid permeability in y direction vector [mDarcy]                      
    reservoir.permz  = np.copy(reservoir.permx)           #fluid permeability in y direction vector [mDarcy]   
    reservoir.Dref   = 7495              #reference depth [feet]
    reservoir.alpha  = 0.0*np.pi/6.0     #dip angle [in radians]

    # Well parameters
    well.x   = [[3637.5],[3787.5],[2500.0],[2575.0],[2650.0],[2200.0],[2275.0],[2350.0],[1087.5],[412.5]]
    well.y   = [[5550.0],[3550.0],[4350.0],[4350.0],[4350.0],[2650.0],[2650.0],[2650.0],[1050.0],[3050.0]]    
    well.type= [[1],[1],[2],[2],[2],[2],[2],[2],[2],[2]] # 1 for rate, 2 for BHP
    well.constraint = [[-120*5.61],[-90*5.61], [502.5], [502.5],[502.5], [502.5],[502.5], [502.5], [502.5], [502.5]]  # rate = scf/day (+ for injector); BHP = psi (always +ve)
    well.rw  = [[0.25],[0.25],[0.25],[0.25],[0.25],[0.25],[0.25],[0.25],[0.25],[0.25]]             # well radius, ft
    well.skin= [[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]] # well skin friction factor, dimensionless
    well.direction = [['v'],['v'],['hx'],['hx'],['hx'],['hx'],['hx'],['hx'],['v'],['v']]   # direction of the well: vertical v or horizontal hx or hy

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
    numerical.D = np.reshape(depth, (numerical.N,1))

    BC.type  = [['Neumann'],['Neumann'],['Neumann'],['Neumann']] #type of BC: left, right, bottom, top
    BC.value = [[0],[0],[0],[0]]    #value of the boundary condition: psi or ft^3/day
    #BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)],[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    #IC.P = 1000.0*np.ones((numerical.N, 1)) #initial condition: pressure in the domain
    IC.P = 4500.0 + (fluid.rho/144.0)*(numerical.D - reservoir.Dref)
depth    =-np.loadtxt("PJ1-Depth.txt")
porosity = np.loadtxt("PJ1-Porosity.txt")
permx    = np.loadtxt("PJ1-Permeability.txt") 

inputfile(fluid,reservoir,numerical,BC,IC,well)

