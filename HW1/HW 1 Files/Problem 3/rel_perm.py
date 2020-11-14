"""
reservoir simulation assignment 1
Capillary pressure and relative permeability: Relative permeability
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/10/2020
"""
import numpy as np

def rel_perm(petro,Sw):
    """s: saturation, s_wp: wetting percolation threshold, s_nwp: non-wetting percolation threshold, mu: viscosity, k_0: relative permeability threshold, n: power law parameter"""
    
    S =(Sw - petro.Swr)/(1.0 - petro.Sor - petro.Swr) #Normalized saturation
    krw = petro.krw0 * S**petro.nw                    #Corey-Brooks model
    kro = petro.kro0 *(1.0-S)**petro.no               #Corey-Brooks model

    #implementing percolation threshold conditions
    #for water phase
    if(Sw <= petro.Swr):
        krw = np.nan
        kro = np.nan
    #for oil phase
    if(Sw >= 1.0-petro.Sor):
        krw = np.nan 
        kro = np.nan       

    return np.array([krw, kro], dtype=np.float64)