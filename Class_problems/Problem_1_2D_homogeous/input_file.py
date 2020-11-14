"""
reservoir simulation assignment 3
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
    
    numerical.dt     = 0.01  #time step (days)
    numerical.tfinal = 0.2   #final time [days]
    numerical.Nx     = 100    #number of grid blocks in x-direction
    numerical.Ny     = 100     #number of grid blocks in y-direction    
    numerical.N  = numerical.Nx * numerical.Ny #Total number of grid blocks
    numerical.theta  =  0.0      #type of method theta = 1-Explicit, 0-Implicit, 0.5-CN
    
    reservoir.L  = 10000.0  #length of the reservoir [feet] 
    reservoir.h  = 20.0     #height of the reservoir [feet] 
    reservoir.W  = 10000.0  #width of the reservoir [feet]    
    reservoir.phi= 0.2*np.ones((numerical.N, 1))   #porosity of the reservior vector
    reservoir.k  = 50.0*np.ones((numerical.N, 1))  #fluid permeability vector [mDarcy] 
    reservoir.Dref=1000/(fluid.rho/144.0)      #reference depth [feet]
    reservoir.alpha=0.0*np.pi/6.0              #dip angle [in radians]
    
    #Defining numerical parameters for discretized solution
    numerical.dx = (reservoir.L/numerical.Nx)*np.ones((numerical.Nx, 1)) #block thickness in x vector
    numerical.dy = (reservoir.W/numerical.Ny)*np.ones((numerical.Ny, 1)) #block thickness in y vector

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
    
    BC.type  = [['Dirichlet'],['Neumann'],['Dirichlet'],['Neumann']] #type of BC: left, right, bottom, top
    BC.value = [[2000],[0],[2000],[0]]    #value of the boundary condition: psi or ft^3/day
    BC.depth = [[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)],[reservoir.Dref],[reservoir.Dref-reservoir.L*np.sin(reservoir.alpha)]]
    
    IC.P = 1000.0*np.ones((numerical.N, 1)) #initial condition: pressure in the domain
    #IC.P = 1000.0 + (fluid.rho/144.0)*(numerical.D - reservoir.Dref)