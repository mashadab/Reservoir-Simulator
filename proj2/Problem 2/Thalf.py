"""
reservoir simulation project 2 (2020): Problem 1
2D reservoir simulation: Interblock transmissibility
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

from rel_perm import rel_perm

#fluid, reservoir and simulation parameters   
def Thalf(i,j,direction,fluid,reservoir,petro,numerical,P,Pc,Sw):

    if direction == 'x':  
        kAd = (2*reservoir.permx[i,0]*numerical.dy[i,0]*reservoir.h*reservoir.permx[j,0]*numerical.dy[j,0]*reservoir.h)/ \
            (reservoir.permx[i,0]*numerical.dy[i,0]*reservoir.h*numerical.dx[j,0] + \
             reservoir.permx[j,0]*numerical.dy[j,0]*reservoir.h*numerical.dx[i,0])
    elif direction == 'y':  
        kAd = (2*reservoir.permy[i,0]*numerical.dx[i,0]*reservoir.h*reservoir.permy[j,0]*numerical.dx[j,0]*reservoir.h)/ \
            (reservoir.permy[i,0]*numerical.dx[i,0]*reservoir.h*numerical.dy[j,0] + \
             reservoir.permy[j,0]*numerical.dx[j,0]*reservoir.h*numerical.dy[i,0])

    #Calculating the potential
    POT_i  = P[i,0]  - ( fluid.rhoo[i,0] / (144.0*fluid.Bo[i,0]) ) * numerical.D[i,0] 
    POT_j  = P[j,0]  - ( fluid.rhoo[j,0] / (144.0*fluid.Bo[j,0]) ) * numerical.D[j,0] 

    if POT_i >= POT_j:
        krw,kro = rel_perm(petro,Sw[i,0])
        Bw  = fluid.Bw[i,0]
        Bo  = fluid.Bo[i,0]    
        muw = fluid.muw[i,0]    
        muo = fluid.muo[i,0]    

    else:
        krw,kro = rel_perm(petro,Sw[j,0])
        Bw  = fluid.Bw[j,0]
        Bo  = fluid.Bo[j,0]    
        muw = fluid.muw[j,0]    
        muo = fluid.muo[j,0]      
    
    fluidhalfw = krw/(muw * Bw)
    fluidhalfo = kro/(muo * Bo)    
    
    Tw = fluidhalfw*kAd
    To = fluidhalfo*kAd
    
    return Tw, To;