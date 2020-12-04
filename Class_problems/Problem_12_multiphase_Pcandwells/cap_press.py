"""
reservoir simulation project 2 (2020): Problem 1
2D Multiphase reservoir simulation: 
Capillary pressure and relative permeability: Capillary pressure
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 12/04/2020
"""
import numpy as np

def cap_press(petro,Sw):
     
    S =(Sw - petro.Swr)/(1.0 - petro.Sor - petro.Swr) #Normalized saturation {Eq. 1.29}
    Se= (Sw - petro.Swr)/(1.0 - petro.Swr)            #{Eq. 1.28a}
    
    #Corey-Brooks model
    Pcd = petro.Pe * Se **(-1.0 / petro.lamda)        #Drainage capillary pressure (used for initialization) {Eq. 1.28a}  

    Pci = petro.Pe * (S **(-1.0 / petro.lamda) - 1.0) #Imbibition capillary pressure (used for initialization) {Eq. 1.28b}
    
    Se2 = Se + 0.001
    Pcd2 = petro.Pe * Se2 **(-1.0 / petro.lamda)
    Pcprime = (Pcd2 - Pcd)/0.001
    
    return Pcd,Pcprime;