"""
reservoir simulation project 2 (2020): Problem 1
2D Multiphase reservoir simulation: Input file
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 11/12/2020
"""

import numpy as np

class fluid:
    def __init__(self):
        self.mu = []
class reservoir:
    def __init__(self):
        self.dt = [] 

#fluid, reservoir parameters   
def fluid_properties(reservoir,fluid,P,Pw):
    fluid.Bw = 1.0*(1.0-fluid.cw*(Pw-1000.0))
    fluid.Bo = fluid.B_BP*(1.0-fluid.co*(P-fluid.BP))
    fluid.rhogsc = (28.9586 * fluid.sg / 379.4)
    fluid.rhoo = (fluid.rhoosc + fluid.Rs*fluid.rhogsc/5.615)/fluid.Bo
    fluid.muo=5 + 0.001*(P-fluid.BP)
    