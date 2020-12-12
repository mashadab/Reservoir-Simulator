"""
reservoir simulation project 2 (2020): Problem 1
2D Multiphase reservoir simulation: 
Capillary pressure and relative permeability: Capillary pressure
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 12/04/2020
"""
import numpy as np

def cap_press(petro,Sw,Sw_hyst):

    S      = (Sw - petro.Swr)/(1.0 - petro.Sor - petro.Swr) #Normalized saturation {Eq. 1.29}
    Se     = (Sw - petro.Swr)/(1.0 - petro.Swr)             #{Eq. 1.28a}
    
    #Corey-Brooks model
    Pcd    = petro.Pe * (Se) **(-1.0 / petro.lamda)        #Drainage capillary pressure (used for initialization) {Eq. 1.28a}  
    Pci    = petro.Pe * ((S) **(-1.0 / petro.lamda) - 1.0) #Imbibition capillary pressure (used for initialization) {Eq. 1.28b}
    
    #Capillary pressure scanning curve
    epspc  = 0.1
    Sw_max = 1.0 - petro.Sor
    f      = ((Sw_max - Sw_hyst + epspc)/(Sw_max - Sw_hyst)) * ((Sw - Sw_hyst) / (Sw - Sw_hyst + epspc))
    #f = 1.0
    Pc     = f * Pci + (1.0 - f) * Pcd

    #Calculate derivative
    S2     = (Sw + 0.001 - petro.Swr) / (1.0 - petro.Swr - petro.Sor)
    Se2    = (Sw + 0.001 - petro.Swr) / (1.0 - petro.Swr)
    
    Pcd2   = petro.Pe * Se2 **(-1.0 / petro.lamda)
    Pci2   = petro.Pe * (S2 **(-1.0 / petro.lamda) - 1.0 )
    f2     = ((Sw_max - Sw_hyst + epspc) / (Sw_max - Sw_hyst)) * ((Sw + 0.001 - Sw_hyst) / (Sw + 0.001 - Sw_hyst + epspc))
    #f2 = 1.0
    Pc2    = f2 * Pci + (1.0 - f2) *  Pcd
    Pcprime = (Pc2 - Pc)/0.001
    
    return Pc,Pcprime;