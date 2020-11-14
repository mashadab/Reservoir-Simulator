"""
reservoir simulation assignment 1
Initializing a reservoir
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/3/2020
"""

import numpy as np
    
def build_grid(Grid):

    # This function computes takes in minimal definition of the computational
    # domain and grid and computes all containing all pertinent information 
    # about the grid. 
    # Input:
    # Grid.Lx = length of the domain
    # Grid.dx = cell width
    # Grid.xc = vector of cell center locations
    # Grid.xf = vector of cell face locations
    # Grid.Nfx = number of fluxes in x-direction
    # Grid.dof_xmin = degrees of fredom corrsponding to the cells along the x-min boundary
    # Grid.dof_xmax = degrees of fredom corrsponding to the cells along the x-max boundary
    # Grid.dof_ymin = degrees of fredom corrsponding to the cells along the y-min boundary
    # Grid.dof_ymax = degrees of fredom corrsponding to the cells along the y-max boundary
    
    # Grid.dof_f_xmin = degrees of fredom corrsponding to the faces at the x-min boundary
    # Grid.dof_f_xmax = degrees of fredom corrsponding to the faces at the x-max boundary
    # Grid.dof_f_ymin = degrees of fredom corrsponding to the faces at the y-min boundary
    # Grid.dof_f_ymax = degrees of fredom corrsponding to the faces at the y-max boundary
    # Example call: 
    # >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10; 
    # >> Grid = build_grid(Grid);
    
    # Set up the geometry
    #In x-direction
    if not hasattr(Grid,'xmin'):
        Grid.xmin = 0
        print("Grid.xmin is not defined and has been set to 0.\n")
    if not hasattr(Grid,'xmax'):
        Grid.xmax = 1 
        print("Grid.xmax is not defined and has been set to 1.\n")
    if not hasattr(Grid,'Nx'): 
        Grid.Nx   = 1
        print("Grid.Nx is not defined and has been set to 1.\n")
        
    Grid.Lx = Grid.xmax-Grid.xmin    # domain length in x
    Grid.dx = Grid.Lx/Grid.Nx        # dx of the gridblocks
    
    #In y-direction
    if not hasattr(Grid,'ymin'):
        Grid.ymin = 0
        print("Grid.ymin is not defined and has been set to 0.\n")
    if not hasattr(Grid,'ymax'):
        Grid.ymax = 1 
        print("Grid.ymax is not defined and has been set to 1.\n")
    if not hasattr(Grid,'Ny'): 
        Grid.Ny   = 1
        print("Grid.Ny is not defined and has been set to 1.\n")
        
    Grid.Ly = Grid.ymax-Grid.ymin    # domain length in y
    Grid.dy = Grid.Ly/Grid.Ny        # dy of the gridblocks
    
    # Number for fluxes
    Grid.Nfx =  (Grid.Nx+1)*Grid.Ny
    Grid.Nfy =  Grid.Nx*(Grid.Ny+1)
    Grid.Nf  =  Grid.Nfx + Grid.Nfy  
    
    # x coords of the corners of the domain
    Grid.xdom = [Grid.xmin,Grid.xmin,Grid.xmax,Grid.xmax]
    Grid.ydom = [Grid.ymin,Grid.ymax,Grid.ymin,Grid.ymax]
    
    #Set up mesh for plotting
    #xcoords of the cell centers    
    #Grid.xc = [Grid.xmin+Grid.dx/2:Grid.dx:Grid.xmax-Grid.dx/2] # x-coords of gridblock centers
    #temp = [xmin+xi*Grid.dx + Grid.dx/2 for xi in range(Grid.Nx)]
    Grid.xc = np.transpose(np.linspace(Grid.xmin+Grid.dx/2, Grid.xmax-Grid.dx/2, Grid.Nx))
    Grid.yc = np.transpose(np.linspace(Grid.ymin+Grid.dy/2, Grid.ymax-Grid.dy/2, Grid.Ny))

    Grid.xf = np.transpose(np.linspace(Grid.xmin, Grid.xmax, Grid.Nx+1)) # x-coords of gridblock faces
    Grid.yf = np.transpose(np.linspace(Grid.ymin, Grid.ymax, Grid.Ny+1)) # y-coords of gridblock faces

    # Set up dof vectors
    Grid.N = Grid.Nx*Grid.Ny                                             # total number of gridblocks
    Grid.dof   = np.transpose([xi+1 for xi in range(Grid.N)])            # cell centered degree of freedom/gridblock number
    Grid.dof_f = np.transpose([xi+1 for xi in range(Grid.Nf)])           # face degree of freedom/face number

    # Boundary dof's
    # Boundary cells
    # make more efficient by avoidng DOF
    DOF = np.transpose(np.reshape(Grid.dof,(Grid.Nx,Grid.Ny)))

    Grid.dof_xmin = DOF[:,0]
    Grid.dof_xmax = DOF[:,Grid.Nx-1]
    Grid.dof_ymin = np.transpose(DOF[0,:])
    Grid.dof_ymax = np.transpose(DOF[Grid.Ny-1,:])

    # Boundary faces
    DOFx = np.transpose(np.array([list(range(1,Grid.Nfx+1,1))]).reshape((Grid.Nx+1,Grid.Ny)))
    Grid.dof_f_xmin = DOFx[:,0]
    Grid.dof_f_xmax = DOFx[:,Grid.Nx+1-1]

    #Grid.dof_f_xmin = Grid.dof_xmin
    #Grid.dof_f_xmax = np.transpose(np.array([list(range(Grid.Nfx-Grid.Ny+1,Grid.Nfx+1,1))]))

    DOFy = np.transpose(np.reshape(Grid.Nfx + np.array([list(range(1,Grid.Nfy+1,1))]),(Grid.Nx,Grid.Ny+1)))
    Grid.dof_f_ymin = np.transpose(DOFy[0,:])
    Grid.dof_f_ymax = np.transpose(DOFy[Grid.Ny+1-1,:])

    Grid.A  = np.concatenate([np.ones((Grid.Nfx,1))*Grid.dy,np.ones((Grid.Nfy,1))*Grid.dx,[[Grid.dx*Grid.dy]]], axis=0 )
    Grid.V  = np.ones((Grid.N,1))*Grid.dx*Grid.dy
    
    return Grid;