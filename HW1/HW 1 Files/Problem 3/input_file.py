"""
reservoir simulation assignment 1
Capillary pressure and relative permeability: Input files
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/10/2020
"""

#Multiphase / Relativity Permeability values
        
def inputfile(petro):
    petro.Swr  = 0.2    #residual water saturation
    petro.Sor  = 0.4    #residual oil saturation
    petro.nw   = 2.0    #Corey-Brooks exponent (water)
    petro.no   = 2.0    #Corey-Brooks exponent (oil)
    petro.krw0 = 0.3    #Corey-Brooks endpoint (water)
    petro.kro0 = 0.8    #Corey-Brooks endpoint (oil)
    
    #Capillary Pressure
    petro.lamda= 2.0     #fitting parameter for Corey-Brooks model
    petro.Pe   = 3.5     #capillary entry pressure [psi]
