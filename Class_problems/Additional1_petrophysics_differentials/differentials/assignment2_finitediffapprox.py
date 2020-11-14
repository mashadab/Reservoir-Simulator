"""
reservoir simulation assignment 2
Pressure values and derivative approximation
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/22/2020
"""

#Importing required libraries
from IPython import get_ipython
get_ipython().magic('reset -sf')#for clearing everything
get_ipython().run_line_magic('matplotlib', 'qt') #for plotting in separate window
# Import curve fitting package from scipy
from scipy.optimize import curve_fit

import numpy as np              #import numpy library
import matplotlib.pyplot as plt #library for plotting
plt.rcParams.update({'font.size': 22})

# Function to calculate the power-law with constants a and b
def power_law(x, a, b):
    return a*np.power(x, b)

#Making grid class
class grid:
    def __init__(self):
        self.xmin = []
        self.xmax = []
        self.Nx = []
        self.N = []

#Input parameters
pi        = 1000        # initial reservior pressure [psi]   
qw        = 100         # well flow rate [ft^3/day]   
k         = 50          # permeability [Darcy]
h         = 20          # thickness of the reservior [feet]
alpha     = 10.0        # diffusivity contant [ft^2/day]   
rw        = 0.25        # well radius [feet]     
mu        = 1.0         # fluid viscosity [centipoise]   
conv_fac  = 1.0/(6.33E-3) #conversion factor for mD per Darcy to ft^2/psi-day
tf        = 1           # final time [day]
dt_a      = 0.01        # time step for approximation part a [day]
dr_b      = 0.001       # radial step for approximation part b [feet]
dr_c      = 0.001       # radial step for approximation part c [feet]

#Calculate the constant  beta
beta      = qw * mu / (4 * np.pi* k * h) *conv_fac 

#Pressure as anonymous function
P = lambda r,t: pi - beta*(np.log(r**2.0/(alpha*t)) + 0.57721)

#Evaluating the analytical values
dPbydt_analy = lambda t: beta / t
dPbydr_analy = lambda r: -2.0 * beta/r
d2Pbydr2_analy = lambda r: 2.0 * beta/r**2.0

#PART B =========================================
#subpart (a)
dPbydt_frwd  = (P(rw,tf+dt_a)-P(rw,tf))/dt_a #forward difference
dPbydt_bkwrd = (P(rw,tf)-P(rw,tf-dt_a))/dt_a #backward difference
dPbydt_cntr  = (P(rw,tf+dt_a)-P(rw,tf-dt_a))/(2*dt_a) #centered difference

#subpart (b)
dPbydr_frwd  = (P(rw+dr_b,tf)-P(rw,tf))/dr_b #forward difference
dPbydr_bkwrd = (P(rw,tf)-P(rw-dr_b,tf))/dr_b #backward difference
dPbydr_cntr  = (P(rw+dr_b,tf)-P(rw-dr_b,tf))/(2*dr_b) #centered difference

#subpart (c)
d2Pbydr2_cntr  = ( P(rw+dr_c,tf) -2.0 * P(rw,tf) + P(rw-dr_c,tf) ) / dr_c**2 #centered difference


#PART C =============================================
dr = np.logspace(-3,-1,100)
dt = np.logspace(-4,-2,100)

#subpart (a)
VdPbydt_frwd  = (P(rw,tf+dt)-P(rw,tf))/dt #forward difference V stands for vectorized
VdPbydt_bkwrd = (P(rw,tf)-P(rw,tf-dt))/dt #backward difference
VdPbydt_cntr  = (P(rw,tf+dt)-P(rw,tf-dt))/(2*dt) #centered difference
#Absolute errors
EdPbydt_frwd  =  np.abs(VdPbydt_frwd - dPbydt_analy(tf))  #errors in forward difference E stands for vectorized errors
EdPbydt_bkwrd = np.abs(VdPbydt_bkwrd - dPbydt_analy(tf))  #errors in backward difference E stands for vectorized errors
EdPbydt_cntr  =  np.abs(VdPbydt_cntr - dPbydt_analy(tf))  #errors in centered difference E stands for vectorized errors
#Power law fit to find the order of accuracy
ZdPbydt_frwd, cov = curve_fit(f=power_law, xdata=dt, ydata=EdPbydt_frwd, bounds=(-np.inf, np.inf))
ZdPbydt_bkwrd, cov = curve_fit(f=power_law, xdata=dt, ydata=EdPbydt_bkwrd, bounds=(-np.inf, np.inf))
ZdPbydt_cntr, cov = curve_fit(f=power_law, xdata=dt, ydata=EdPbydt_cntr, bounds=(-np.inf, np.inf))

#subpart (b)
VdPbydr_frwd  = (P(rw+dr,tf)-P(rw,tf))/dr #forward difference
VdPbydr_bkwrd = (P(rw,tf)-P(rw-dr,tf))/dr #backward difference
VdPbydr_cntr  = (P(rw+dr,tf)-P(rw-dr,tf))/(2*dr) #centered difference
#Absolute errors
EdPbydr_frwd  =  np.abs(VdPbydr_frwd - dPbydr_analy(rw))  #errors in forward difference E stands for vectorized errors
EdPbydr_bkwrd  = np.abs(VdPbydr_bkwrd - dPbydr_analy(rw)) #errors in backward difference E stands for vectorized errors
EdPbydr_cntr  =  np.abs(VdPbydr_cntr - dPbydr_analy(rw))  #errors in centered difference E stands for vectorized errors

#Power law fit to find the order of accuracy
ZdPbydr_frwd, cov = curve_fit(f=power_law, xdata=dr, ydata=EdPbydr_frwd, bounds=(-np.inf, np.inf))
ZdPbydr_bkwrd, cov = curve_fit(f=power_law, xdata=dr, ydata=EdPbydr_bkwrd, bounds=(-np.inf, np.inf))
ZdPbydr_cntr, cov = curve_fit(f=power_law, xdata=dr, ydata=EdPbydr_cntr, bounds=(-np.inf, np.inf))

#subpart (c)
Vd2Pbydr2_cntr  = (P(rw+dr,tf) -2.0 * P(rw,tf) + P(rw-dr,tf) ) / dr**2 #centered difference
#Absolute errors
Ed2Pbydr2_cntr  =  np.abs(Vd2Pbydr2_cntr - d2Pbydr2_analy(rw))  #errors in centered difference E stands for vectorized errors
#Polyfit to find the order of accuracy
Zd2Pbydr2_cntr, cov = curve_fit(f=power_law, xdata=dr, ydata=Ed2Pbydr2_cntr, bounds=(-np.inf, np.inf))

#PLOTTING
#Higher order derivative approximations ===========================================================
#(a) Plotting the forward, backward and center difference for dP/dt
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(dt,VdPbydt_frwd,'ro',label=r'$Forward$')
plot = plt.plot(dt,VdPbydt_bkwrd,'bd',label=r'$Backward$')
plot = plt.plot(dt,VdPbydt_cntr,'gs',label=r'$Centred$')
plt.hlines(dPbydt_analy(tf),dt[1],dt[99],label=r'$Analytical$')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$dP/dt$ [psi/day]')
plt.xlabel(r'$\Delta t$ [day]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'dpbydt.png',bbox_inches='tight', dpi = 600)

#(b) Plotting the forward, backward and center difference for dP/dr
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(dr,VdPbydr_frwd,'ro',label=r'$Forward$')
plot = plt.plot(dr,VdPbydr_bkwrd,'bd',label=r'$Backward$')
plot = plt.plot(dr,VdPbydr_cntr,'gs',label=r'$Centred$')
plt.hlines(dPbydr_analy(rw),dr[1],dr[99],label=r'$Analytical$')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$dP/dr$ [psi/feet]')
plt.xlabel(r'$\Delta r$ [feet]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'dpbydr.png',bbox_inches='tight', dpi = 600)

#(b) Plotting the center difference for d^2P/dr^2
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.plot(dr,Vd2Pbydr2_cntr,'ro',label=r'$Centred$')
plt.hlines(d2Pbydr2_analy(rw),dr[1],dr[99],label=r'$Analytical$')
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$d^2P/dr^2$ [psi/feet^2]')
plt.xlabel(r'$\Delta r$ [feet]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'd2pbydr2.png',bbox_inches='tight', dpi = 600)

#ABSOLUTE ERRORS===========================================================
#(a) Plotting the absolute errors in forward, backward and center difference for dP/dt
ax = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.loglog(dt,EdPbydt_frwd,'ro',label=r'$Forward$')
plot = plt.loglog(dt,power_law(dt, *ZdPbydt_frwd),'r--',label=r'$fit: O( \Delta t^{%5.3f})$' %ZdPbydt_frwd[1] )
plot = plt.loglog(dt,EdPbydt_bkwrd,'bd',label=r'$Backward$')
plot = plt.loglog(dt,power_law(dt, *ZdPbydt_bkwrd),'b--',label=r'$fit: O(\Delta t^{%5.3f})$' %ZdPbydt_bkwrd[1] )
plot = plt.loglog(dt,EdPbydt_cntr,'gs',label=r'$Centred$')
plot = plt.loglog(dt,power_law(dt, *ZdPbydt_cntr),'g--',label=r'$fit: O(\Delta t^{%5.3f})$' %ZdPbydt_cntr[1] )
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$|dP/dt-dP/dt_{analytical}|$ [psi/day]')
plt.xlabel(r'$\Delta t$ [day]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'AbsErr_dpbydt.png',bbox_inches='tight', dpi = 600)

#(b) Plotting the absolute errors in forward, backward and center difference for dP/dr
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.loglog(dr,EdPbydr_frwd,'ro',label=r'$Forward$')
plot = plt.loglog(dr,power_law(dr, *ZdPbydr_frwd),'r--',label=r'$fit: O( \Delta r^{%5.3f})$' %ZdPbydr_frwd[1] )
plot = plt.loglog(dr,EdPbydr_bkwrd,'bd',label=r'$Backward$')
plot = plt.loglog(dr,power_law(dr, *ZdPbydr_bkwrd),'b--',label=r'$fit: O( \Delta r^{%5.3f})$' %ZdPbydr_bkwrd[1] )
plot = plt.loglog(dr,EdPbydr_cntr,'gs',label=r'$Centred$')
plot = plt.loglog(dr,power_law(dr, *ZdPbydr_cntr),'g--',label=r'$fit: O(\Delta r^{%5.3f})$' %ZdPbydr_cntr[1] )
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.xlabel(r'$\Delta r$ [feet]')
plt.ylabel(r'$|dP/dr-dP/dr_{analytical}|$ [psi/feet]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'AbsErr_dpbydr.png',bbox_inches='tight', dpi = 600)

#(c) Plotting the absolute errors in center difference for d^2P/dr^2
fig = plt.figure(figsize=(15,7.5) , dpi=100)
plot = plt.loglog(dr,Ed2Pbydr2_cntr,'gs',label=r'$Centred$')
plot = plt.loglog(dr,power_law(dr, *Zd2Pbydr2_cntr),'g--',label=r'$fit: O(\Delta r^{%5.3f})$' %Zd2Pbydr2_cntr[1] )
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.ylabel(r'$|d^2P/dr^2-d^2P/dr^2_{analytical}|$ [psi/feet^2]')
plt.xlabel(r'$\Delta r$ [feet]')
legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f'AbsErr_d2pbydr2.png',bbox_inches='tight', dpi = 600)
