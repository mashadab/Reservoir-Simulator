"""
reservoir simulation assignment 1
Capillary pressure and relative permeability: Petrophysics main function
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/10/2020
"""

#import inbuilt libraries
import numpy as np

#importing personal libraries
from input_file import inputfile
from rel_perm import rel_perm
from petroplots import petroplots
from cap_press import cap_press

#making petrophysics class
class petro:
    def __init__(self):
        self.Sor = []

inputfile(petro)     #uploading petrophysical properties           

Sw = np.transpose([np.linspace(0,1,10000)]) #initializing water saturation as a column vector
#Sw_actual = np.transpose([np.linspace(petro.Swr,1-petro.Sor,10000)]) #initializing water saturation

krw = np.zeros((len(Sw),1)) #allocating column vector for relative permeability (water)
kro = np.zeros((len(Sw),1)) #allocating column vector for relative permeability (oil)
Pci = np.zeros((len(Sw),1)) #allocating column vector for imbibition capillary pressure [psi]
Pcd = np.zeros((len(Sw),1)) #allocating column vector for drainage capillary pressure  [psi]

#evaluating relative permeabilities using Corey-Brooks model
for j in range(0, len(Sw)):
    [krw[j], kro[j]] = rel_perm(petro,Sw[j,0])   
    [Pci[j], Pcd[j]] = cap_press(petro,Sw[j,0]) 

#plotting
petroplots(Sw,krw,kro,Pcd,Pci)