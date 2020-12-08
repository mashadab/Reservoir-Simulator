"""
reservoir simulation project 2 (2020): Problem 1
2D Multiphase reservoir simulation: 
Capillary pressure and relative permeability: Rel perm
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 11/12/2020
"""

import numpy as np

def rel_perm(petro,Sw):
    """s: saturation, s_wp: wetting percolation threshold, s_nwp: non-wetting percolation threshold, mu: viscosity, k_0: relative permeability threshold, n: power law parameter"""
    
    S =(Sw - petro.Swr)/(1.0 - petro.Sor - petro.Swr) #Normalized saturation
    
    krw = (petro.krw0 * S**petro.nw)                    #Corey-Brooks model
    kro = (petro.kro0 *(1.0-S)**petro.no)               #Corey-Brooks model

    return np.array([krw, kro], dtype=np.float64)