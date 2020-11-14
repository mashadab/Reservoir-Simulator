"""
reservoir simulation assignment 1
Capillary pressure and relative permeability
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/10/2020
"""

import numpy as np
import matplotlib.pyplot as plt

def mobility(s,s_wp,s_nwp,mu_w,mu_nw,k_w0,k_nw0,n_w,n_nw): 
    """s: saturation, s_wp: wetting percolation threshold, s_nwp: non-wetting percolation threshold, mu: viscosity, k_0: relative permeability threshold, n: power law parameter"""
    
    s_star=(s-s_wp)/(1.0-s_nwp-s_wp)
    kr_w = k_w0*s_star**n_w
    kr_nw= k_nw0*(1.0-s_star)**n_nw

    """percolation threshold conditions"""
    
    if(s<s_wp):
        kr_w = 0.0
        kr_nw= k_nw0

    if(s>1-s_nwp):
        kr_w = k_w0 
        kr_nw= 0.0       
    
    lambda_w = kr_w/mu_w
    lambda_nw= kr_nw/mu_nw
    lambda_t = lambda_w+lambda_nw
    
    f_w      = lambda_w/lambda_t
    f_nw     = lambda_nw/lambda_t
    return np.array([s, kr_w, kr_nw, lambda_w, lambda_nw,lambda_t, f_w, f_nw], dtype=np.float64)


'''
#main program
N     = 1001 # of points
s_wp  = 0.0  #percolation threshold: wetting phase
s_nwp = 0.0  #percolation threshold: non-wetting phase
mu_w  = 1.0  #dynamic viscosity: wetting phase    
mu_nw = 1.0  #dynamic viscosity: non-wetting phase   
k_w0  = 1.0  #relative permeability threshold: wetting phase   
k_nw0 = 1.0  #relative permeability threshold: non-wetting phase   
n_w   = 1.0  #power law coefficient: wetting phase  
n_nw  = 1.0  #power law coefficient: non-wetting phase   

s=np.linspace(start = 0, stop = 1, num = N)

mobility_data = np.array([], dtype=np.float64)
mobility_data = mobility_data.reshape(0,8)
for i in range(0, N):
    addition = mobility(s[i],s_wp,s_nwp,mu_w,mu_nw,k_w0,k_nw0,n_w,n_nw)
    mobility_data = np.vstack([mobility_data, addition])
    #pass

#Plotting
fig, ax= plt.subplots()
#fig.add_axes([0,0,1.0,1.0])
ax.plot(mobility_data[:,0],mobility_data[:,1],'r-',label=r'$k_{rw}$')
ax.plot(mobility_data[:,0],mobility_data[:,2],'b-',label=r'$k_{rnw}$')
legend = ax.legend(loc='upper center', shadow=False, fontsize='x-large')
ax.set_xlabel(r'$s_w$')
ax.set_ylabel(r'$k_r$')
plt.savefig('mobility.png',bbox_inches='tight', dpi = 600)
'''