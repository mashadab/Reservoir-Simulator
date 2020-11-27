"""
reservoir simulation assignment 4
1D reservoir simulation Q1.4: Input files
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, identity

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
class well:
    def __init__(self):
        self.dt = [] 

#fluid, reservoir and simulation parameters   
def inputfile():

  # Numerical simulation properties 
    numerical.Nx = 30       #number of grid blocks in x-direction
    numerical.Ny = 30       #number of grid blocks in y-direction 
    numerical.N  = numerical.Nx * numerical.Ny #Total number of grid blocks
    

  # Reservoir parameters    
    reservoir.L  = 1200.0  #length of the reservoir [feet] 
    reservoir.h  = 20.0    #height of the reservoir [feet] 
    reservoir.W  = 1200.0  #width of the reservoir [feet]        
    reservoir.phi   = 0.20*np.ones((numerical.N, 1))    #porosity of the reservior vector [mDarcy] 
    reservoir.perm  = 100.0* np.ones((numerical.N, 1))  #fluid permeability in x direction vector [mDarcy] 
      
    
    #Fluid properties    
    fluid.mu  = 1.0*np.ones((numerical.N, 1))   #fluid viscosity [centipoise]  
    fluid.Bw  = 1.0*np.ones((numerical.N, 1))   #formation volume factor of water [rb/stb]  
    fluid.cf  = 5E-6        #fluid compressibility [1/psi] 
    fluid.rho = 62.4        #density of the fluid [lbm/ft^3]
    fluid.relperm = 1.0*np.ones((numerical.N, 1))  #Relative permeability of the fluid


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
    
    numerical.x  = np.reshape(numerical.Xc, (numerical.N,1))    #building the single X column vector
    numerical.y  = np.reshape(numerical.Yc, (numerical.N,1))    #building the single Y column vector
    
    #Dimensionless properties
    N_Pe             = np.array([100])
    numerical.type   = 'profile'        #Type 'profile' or 'exit'
    
    #Discretized solution parameters
    numerical.pie    = N_Pe/numerical.N #Block Peclet number (u*dx)/(D*phi)
    numerical.cr     = np.array([0.5])  #Dimensionless Courant number
    #numerical.eta    = numerical.cr/numerical.pie #Dimensionless block size
    numerical.eta    = 0.001            #Dimensionless block size
    
    numerical.dxD    = 1/numerical.N    #Dimensionless block size
    numerical.PVI    = 2.0

    numerical.xD      = np.empty((numerical.N,1))
    #position of the block centres
    for i in range(0,numerical.N):
        numerical.xD[i,0] = ((1.0/numerical.N)*(i+1.0)) - (1.0/numerical.N)/2
        
    numerical.method = 'explicit'  #'Implicit' or 'explicit'
    numerical.scheme = 'upwind'    #'upwind' or 'center'

    #numerical.dx     = numerical.L/numerical.N
    #numerical.dt     = 2.1064e3*3/(numerical.N**2 / 25)
    numerical.dt      = 0.1  
    
    #Initial condition
    IC.C = np.zeros((numerical.N, 1))  #Concentration
    
    # Well parameters
    well.x   = [[1.0],[0.99*reservoir.L]]
    well.y   = [[1.0],[0.99*reservoir.W]]    
    well.type= [[1],[1]]                        # 1 for rate, 2 for BHP
    well.constraint = [[5000],[-5000]]  # rate = scf/day (+ for injector); BHP = psi (always +ve)
    well.rw  = [[0.5],[0.5]]         # well radius, ft
    well.skin= [[0],[0]]             # well skin friction factor, dimensionless
    well.direction = [['v'],['v']]   # direction of the well: vertical v or horizontal hx or hy
    #well.block = [[1],[9]]           # direction of the well: vertical v or horizontal hx or hy

    #Making TAO matrix
    TAO = lil_matrix((numerical.N, numerical.N))
    
    #Adding wells
    TAO[0,0] = well.constraint[0][0]*numerical.dt/(numerical.dx[0,0]*numerical.dy[0,0]*reservoir.h*reservoir.phi[0,0])
    TAO[numerical.N-1,numerical.N-1] = well.constraint[1][0]*numerical.dt/(numerical.dx[0,0]*numerical.dy[0,0]*reservoir.h*reservoir.phi[0,0])
     
    IC.P = np.loadtxt("Pss2.txt")
    IC.P = np.transpose(np.array([IC.P]))
    
    well_type = 1
    return numerical, reservoir, fluid, IC, N_Pe, TAO, well, well_type

numerical, reservoir, fluid, IC, N_Pe, TAO, well, well_type = inputfile()

