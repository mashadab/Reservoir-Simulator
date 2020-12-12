"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Plotting
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

#Importing required libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
class fluid:
    def __init__(self):
        self.mu = []
class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = [] 

#fluid, reservoir and simulation parameters   
def initial_plot(reservoir,numerical,P,time,PorSw):
    if PorSw == 'P':
        fig = plt.figure(figsize=(15,7.5) , dpi=100)
        ax1= plt.contourf(numerical.Xc, numerical.Yc,np.reshape(P,(numerical.Ny,numerical.Nx)),100,cmap=cm.coolwarm, antialiased=True)
        # This is the fix for the white lines between contour levels
        for c in ax1.collections:
            c.set_edgecolor("face")
        plt.title('Pressure at %0.1f [days]' %time)
        plt.axis('scaled')
        plt.xlabel(r'$x [feet]$')
        plt.ylabel(r'$y [feet] $')
        clb = plt.colorbar(ax1)
        clb.set_label(r'$Pressure [psi]$', labelpad=-40, y=1.1, rotation=0)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(f'reservoir_pressure{time}days.png',bbox_inches='tight')
        
    if PorSw == 'Sw':
        fig = plt.figure(figsize=(15,7.5) , dpi=100)
        ax1= plt.contourf(numerical.Xc, numerical.Yc,np.reshape(P,(numerical.Ny,numerical.Nx)),100,cmap=cm.coolwarm, antialiased=True)
        # This is the fix for the white lines between contour levels
        for c in ax1.collections:
            c.set_edgecolor("face")
        plt.title('Water saturation at %0.1f [days]' %time)
        plt.axis('scaled')
        plt.xlabel(r'$x [feet]$')
        plt.ylabel(r'$y [feet] $')
        clb = plt.colorbar(ax1)
        clb.set_label(r'$Sw$', labelpad=-40, y=1.1, rotation=0)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(f'Sw{time}days.png',bbox_inches='tight')
