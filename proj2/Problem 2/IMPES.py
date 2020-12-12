"""
reservoir simulation assignment 10
1D reservoir simulation Q10: Main file (Implicit pressure explicit saturation)
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 12/4/2020
"""

#import inbuilt libraries
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, identity
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
from scipy.sparse.construct import eye
import matplotlib.pyplot as plt
import time as timer
from math import floor, ceil
import warnings
warnings.filterwarnings("ignore")

#importing personal libraries
from input_file_2D import inputfile
from myarrays import myarrays
from updatewells import updatewells
from rel_perm import rel_perm
from spdiaginv import spdiaginv
from init_plot import initial_plot
from postprocess import postprocess

#making simulation and IC classes
class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = [] 
class fluid:
    def __init__(self):
        self.dt = [] 
class grid:
    def __init__(self):
        self.xmin = []
class BC:
    def __init__(self):
        self.xmin = []
class IC:
    def __init__(self):
        self.xmin = []
class petro:
    def __init__(self):
        self.xmin = []
class well:
    def __init__(self):
        self.xmin = []

tprogstart= timer.clock()

#loading inputfile
inputfile(fluid,reservoir,petro,numerical,BC,IC,well)
#Implicit pressure and explicit saturation for update
#looping through time
nmax  = ceil(numerical.tfinal / numerical.dt)
t = np.empty((nmax + 1))  #dimensional time
t[0]= 0
t_D = np.empty((nmax + 1))#non dimensional time
t_D[0]= 0
k     = 0
PV    = 0
P     = np.copy(IC.P)
Pw    = np.copy(IC.Pw)
Sw    = np.array(np.copy(IC.Sw))
Sw_hyst=np.zeros((numerical.N,2))
Sw_hyst[:,0]=Sw[:,0]

fw   = np.empty((nmax +1 ))          #fractional flow of wetting phase
fw[0]= 0
P_plot= np.zeros((numerical.N,nmax + 1)) #matrix to save pressure 
P_plot[:,0] = IC.P[:,0] 
Sw_plot= np.zeros((numerical.N, nmax + 1)) #matrix to save pressure 
Sw_plot[:,0]= IC.Sw[:,0] 
well.typetime       = np.kron(np.ones((nmax,1)),np.transpose(well.type))
well.constrainttime = np.kron(np.ones((nmax,1)),np.transpose(well.constraint))
well.fw             = np.zeros((len(well.x),nmax+1))
well.Qwf            = np.zeros((len(well.x),nmax+1))
well.Jwvectime      = np.zeros((len(well.x),nmax+1))
well.Jovectime      = np.zeros((len(well.x),nmax+1))
well.Jvectime       = np.zeros((len(well.x),nmax+1))
well.Q              = np.zeros((len(well.x),nmax+1))

while (t[k] < numerical.tfinal): #non dimensional time marching    
    if(k == nmax/4 or k == nmax/2 or k == nmax*3/4 or k == nmax): print(k, t[k],Sw,P)
    
    P_old = np.copy(P)   #Placeholdering the old array
    Sw_old= np.copy(Sw)   #Placeholdering the old array

    #Switching according to the problem
    ##############################################################################
    #Switching after 500 days 
    if t[k] > numerical.tswitch:
        well.type= [[2],[2],[2],[2],[2],[2],[2],[2],[1],[1]] # 1 for rate, 2 for BHP
        well.constraint = [[502.5],[502.5],[502.5], [502.5], [502.5], [502.5], [502.5], [502.5],[500*5.61], [500*5.61]]  # rate = scf/day (+ for injector); BHP = psi (always +ve)

    elif k > 0:
            #Well 1
        if well.constraint[0][0] > 0: #injector well
            # Well 1 BHP
            if P_plot[well.block[0][0],k] + well.constrainttime[k][0] / (well.Jovectime[0][k]) <= 502.5:
                well.type[0][0] = 2
                well.constraint[0][0] = 502.5
        
        else: #producer well
            # Well 1 BHP
            if P_plot[well.block[0][0],k] + (1 - well.fw[0][k]) * well.constrainttime[k][0] / (well.Jovectime[0][k]) <= 502.5:
                well.type[0][0] = 2
                well.constraint[0][0] = 502.5

            # Well 2
        if well.constraint[1][0] > 0: #injector well
            # Well 2 BHP
            if P_plot[well.block[1][0],k] + well.constrainttime[k][1] / (well.Jovectime[1][k]) <= 502.5:
                well.type[1][0] = 2
                well.constraint[1][0] = 502.5
        
        else: #producer well
            # Well 2 BHP
            if P_plot[well.block[1][0],k] + (1 - well.fw[1][k]) * well.constrainttime[k][1] / (well.Jovectime[1][k]) <= 502.5:
                well.type[1][0] = 2
                well.constraint[1][0] = 502.5
    
    well.typetime[k,:] = np.copy(np.transpose(well.type[:][:]))
    well.constrainttime[k,:] = np.copy(np.transpose(well.constraint[:][:]))
    ###############################################    

    #Calculating the arrays
    Tw, To, T, d11, d12, d21, d22, D, G, Pc, Pw = myarrays(fluid,reservoir,petro,numerical,BC,P,Sw,Sw_hyst)  

    #updating the wells
    well, Qw, Qo, Jw, Jo = updatewells(reservoir,fluid,numerical,petro,P,Sw,well)

    J = -d22 @ ( spdiaginv(d12) @ Jw ) + Jo

    Q = -d22 @ ( spdiaginv(d12) @ Qw ) + Qo + reservoir.Pwf * J @ np.ones((numerical.N,1))   #Pwf = 800 psi

    if numerical.method == 'IMPES':
        IM = T + J + D          #implicit part coefficient in Eq. 3.44   
        EX = D @ P_old + Q + G  #explicit part or RHS of Eq. 3.44
    
        P = np.transpose([spsolve(IM,EX)]) #solving IM*P = EX or Ax=B        
        Sw = Sw + spdiaginv(d12) @ (-Tw @ (P - (fluid.rhow/144.0) * numerical.D - Pc) - d11 @ (P - P_old) + Qw + Jw @ (reservoir.Pwf - P))         #explicit saturation

    Sw[Sw > 1.0] = 1.0
    Sw[Sw < 0.2] = 0.2

    for i in range(0, numerical.N):
        if Sw[i,0] > Sw_old[i,0] and Sw_hyst[i,1] == 0:  # [i,1] is a flag
            Sw_hyst[i,0] = Sw[i,0]
            Sw_hyst[i,1] = 1.0

        elif Sw[i,0] < Sw_old[i,0]:
            Sw_hyst[i,0] = Sw[i,0]          

    k = k+1    
    P_plot[:,k] = P[:,0]  
    Sw_plot[:,k]= np.array(Sw)[:,0]
    t[k]= t[k-1] + numerical.dt
    t_D[k]= well.constraint[0][0]*t[k-1]/(reservoir.L*reservoir.W*reservoir.h*reservoir.phi[0,0])

    well.Qwf[:,k]  = (Qw[well.block,0].toarray())[:,0] + (Qo[well.block,0].toarray())[:,0] + (well.Jwvec[:,0] + well.Jovec[:,0]) * reservoir.Pwf

    for i in range(0,len(well.x)):
        kblock  = well.block[i][0]
        krw,kro = rel_perm(petro,Sw[kblock,0])
        M = (kro*fluid.muw[kblock,0])/(krw*fluid.muo[kblock,0])
        well.fw[i,k] = 1/(1+M)
        well.Jwvectime[i,k] = well.Jwvec[i,0]
        well.Jovectime[i,k] = well.Jovec[i,0]
        well.Jvectime[i,k]  = J[kblock,kblock]
        well.Q[i,k]         = Q[kblock,0]

    #Setting the well constraints for Well 1 and Well 2 back to normal for t < 500 days
    if t[k] < numerical.tswitch:
        well.type[0][0]       =  1           #Well 1
        well.constraint[0][0] = -120*5.61
        well.type[1][0]       =  1           #Well 2
        well.constraint[1][0] = -90*5.61

P_plot[np.argwhere(reservoir.permx < 0.0001)] = np.nan

tprogend= timer.clock()
print('Time elapsed in the program', tprogend - tprogstart)

np.savez(f'Project2_n{numerical.N}_k{k}', P_plot = P_plot, Sw_plot = Sw_plot, Nx = numerical.Nx, Ny = numerical.Ny,fw =fw,t = t, x1 = numerical.x1, y1 = numerical.y1)

#post process
P_plot[np.argwhere(numerical.D==0),:] = np.nan
Sw_plot[np.argwhere(numerical.D ==0),:] = np.nan

'''
#Create the plots
initial_plot(reservoir,numerical,P_plot[:,0],t[0],'P')
initial_plot(reservoir,numerical,P_plot[:,500],t[500],'P')
initial_plot(reservoir,numerical,P_plot[:,750],t[750],'P')
initial_plot(reservoir,numerical,P_plot[:,nmax],t[nmax],'P')

initial_plot(reservoir,numerical,Sw_plot[:,0],t[0],'Sw')
initial_plot(reservoir,numerical,Sw_plot[:,500],t[500],'Sw')
initial_plot(reservoir,numerical,Sw_plot[:,750],t[750],'Sw')
initial_plot(reservoir,numerical,Sw_plot[:,nmax],t[nmax],'Sw')

postprocess(P_plot,numerical,well,t)
'''