"""
reservoir simulation assignment 1
Initializing a reservoir
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/3/2020
"""

#Importing required libraries

import numpy as np              #import numpy library
import matplotlib.pyplot as plt #library for plotting
import build_gridfun2D          #personal library for making grid
plt.rcParams.update({'font.size': 22})
from matplotlib import cm

#Making grid class
class grid:
    def __init__(self):
        self.xmin = []
        self.xmax = []
        self.Nx = []
        self.N = []

#Input parameters
length    = 6000        # length of the reservoir [feet]   
breadth   = 7500        # breadth of the reservoir [feet]   
Nx        = 80          # cells in x-direction
Ny        = 75          # cells in y-direction
thic      = 50          # thickness of the reservoir [feet]   
D_woc     = 7475.45     # depth of water oil contact line [feet]     
P_w_woc   = 4500        # pressure at water oil contact line [psi]   
rho_w_sc  = 62.4        # density of water at standard conditions [lbm/ft^3] 
rho_o_sc  = 53.0        # density of oil at standard conditions [lbm/ft^3]
rho_g_sc  = 0.0458171   # density of gas at standard conditions [lbm/ft^3]
c_w       = 2.87e-6     # compressibility of water [1/psi] 
c_o       = 3e-5        # compressibility of oil [1/psi]  
B_w       = 1.012298811 # formation volume factor of water 
B_o       = 1.04567     # formation volume factor of oil 
Rs        = 90.7388     # solution gas oil ratio [ft^3/bbl]
Pe        = 3.5         # capillary entry pressure at water oil contact line [psi]   
s_wr      = 0.2         # residual water saturation
s_or      = 0.4         # residual oil saturation
lam       = 2           # model parameter for capillary pressure

#Calculate water and oleic phase densities
rho_w     = rho_w_sc / B_w
rho_o     = (rho_o_sc + rho_g_sc * Rs /5.61)/B_o #5.61 ft^3 in a barrel (bbl)
    
#Building grid
grid.xmin = 0.0
grid.xmax = length
grid.Nx   = Nx

grid.ymin = 0.0
grid.ymax = breadth
grid.Ny   = Ny

build_gridfun2D.build_grid(grid)                       #building grid
[Xc,Yc] = np.meshgrid(grid.xc,grid.yc)                 #building the (x,y) matrix

#Importing depth text file
depth = np.loadtxt("Inc_Depth.dat")
depth[depth==0.0] = np.nan                             #removing the zero depth for nice plotting

# Variation of pressure and saturation with depth (ignoring compressibility)
P_w = P_w_woc + rho_w / 144.0 * (depth - D_woc) #144 in^2 in 1 ft^2 
P_o = (P_w_woc + Pe) + rho_o / 144.0 * (depth - D_woc)  
s_w = s_wr + (1-s_wr)*((P_o-P_w)/Pe) **(-lam)  #From Corey-Brooks draining curve
s_w[depth>=D_woc] = 1.0
s_o = 1-s_w

#converting to a column vector
depth_col= np.reshape(np.transpose(depth), (grid.N,-1))  #building the single depth vector
P_o_col  = np.reshape(np.transpose(P_o), (grid.N,-1))    #building the single P_o vector
P_w_col  = np.reshape(np.transpose(P_w), (grid.N,-1))    #building the single P_w vector
s_w_col  = np.reshape(np.transpose(s_w), (grid.N,-1))    #building the single s_w vector
s_o_col  = np.reshape(np.transpose(s_o), (grid.N,-1))    #building the single s_o vector

#Plotting starts here

#(a) Plotting the pressure variation with depth
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(P_w_col,depth_col,'b.',label=r'$P_{water}$ [psi]')
plot = plt.plot(P_o_col,depth_col,'r.',label=r'$P_{oil}$ [psi]')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$Pressure$ [psi]')
plt.ylabel(r'$Depth, D$ [m]')
plt.gca().invert_yaxis()
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'Pressures.png',bbox_inches='tight', dpi = 600)

#(b) Plotting the water saturation variation with depth
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(s_w_col,depth_col,'b.')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$s_{w}$')
plt.ylabel(r'$Depth, D$ [m]')
plt.gca().invert_yaxis()
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'Saturation.png',bbox_inches='tight', dpi = 600)

#(c) Plotting the 2D contour of depth
fig = plt.figure(figsize=(15,7.5) , dpi=100)
ax1= plt.contourf(Xc, Yc,depth,cmap=cm.coolwarm, antialiased=True)
plt.title('Depth of the reservoir')
plt.axis('scaled')
plt.xlabel(r'$x [feet]$')
plt.ylabel(r'$y [feet] $')
clb = plt.colorbar(ax1)
clb.set_label(r'$Depth,D [feet]$', labelpad=-40, y=1.1, rotation=0)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'reservoir_depth.png',bbox_inches='tight', dpi = 600)

#(d) Plotting the 2D contour of water saturation
fig = plt.figure(figsize=(15,7.5) , dpi=100)
ax1= plt.contourf(Xc, Yc,s_w,cmap=cm.coolwarm, antialiased=True)
plt.title('Water saturation in the reservoir')
plt.axis('scaled')
plt.xlabel(r'$ x [feet]$')
plt.ylabel(r'$y [feet] $')
clb = plt.colorbar(ax1)
clb.set_label(r'$s_w$', labelpad=-40, y=1.1, rotation=0)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'reservoir_water_saturation.png',bbox_inches='tight', dpi = 600)

#(e) Plotting the 2D contour of oil pressure
fig = plt.figure(figsize=(15,7.5) , dpi=100)
ax1= plt.contourf(Xc, Yc,P_o,cmap=cm.coolwarm, antialiased=True)
plt.title('Oil pressure in the reservoir')
plt.axis('scaled')
plt.xlabel(r'$x [feet]$')
plt.ylabel(r'$y [feet] $')
clb = plt.colorbar(ax1)
clb.set_label(r'$P_{oil} [psi]$', labelpad=-40, y=1.1, rotation=0)
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'reservoir_oil_pressure.png',bbox_inches='tight', dpi = 600)