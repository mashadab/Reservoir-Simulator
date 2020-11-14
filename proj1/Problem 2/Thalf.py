"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Interblock transmissibility
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

#fluid, reservoir and simulation parameters   
def Thalf(i,j,direction,fluid,reservoir,numerical):
 
    if direction == 'x':  
        kAd = (2*reservoir.permx[i,0]*numerical.dy[i,0]*reservoir.h*reservoir.permx[j,0]*numerical.dy[j,0]*reservoir.h)/ \
            (reservoir.permx[i,0]*numerical.dy[i,0]*reservoir.h*numerical.dx[j,0] + \
             reservoir.permx[j,0]*numerical.dy[j,0]*reservoir.h*numerical.dx[i,0])
    elif direction == 'y':  
        kAd = (2*reservoir.permy[i,0]*numerical.dx[i,0]*reservoir.h*reservoir.permy[j,0]*numerical.dx[j,0]*reservoir.h)/ \
            (reservoir.permy[i,0]*numerical.dx[i,0]*reservoir.h*numerical.dy[j,0] + \
             reservoir.permy[j,0]*numerical.dx[j,0]*reservoir.h*numerical.dy[i,0])
    fluidhalf = 0.5*fluid.relperm[i,0]/(fluid.mu[i,0]*fluid.Bw[i,0]) + 0.5*fluid.relperm[j,0]/(fluid.mu[j,0]*fluid.Bw[j,0])
    T = fluidhalf*kAd
    return T;