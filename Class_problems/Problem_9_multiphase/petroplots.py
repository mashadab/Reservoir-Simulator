"""
reservoir simulation assignment 1
Capillary pressure and relative permeability: Petrophysics plots
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/10/2020
"""
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})

def petroplots(Sw,krw,kro,Pcd,Pci):

    #Plotting
    fig, ax= plt.subplots()
    ax.plot(Sw[:,0],krw[:,0],'b-',label=r'$kr_w$')
    ax.plot(Sw[:,0],kro[:,0],'r-',label=r'$kr_o$')
    ax.legend(loc='best', shadow=False, fontsize='x-large')
    ax.set_xlabel(r'Water saturation, $Sw$')
    ax.set_ylabel(r'Relative permeability, $kr$')
    plt.xlim([0,1])
    plt.savefig('relperm.png',bbox_inches='tight', dpi = 600)
    
    #Plotting
    fig, ax= plt.subplots()
    ax.plot(Sw[:,0],Pcd[:,0],'b-',label=r'$Pc_d$')
    ax.plot(Sw[:,0],Pci[:,0],'r-',label=r'$Pc_i$')
    ax.legend(loc='best', shadow=False, fontsize='x-large')
    ax.set_xlabel(r'Water saturation, $Sw$')
    ax.set_ylabel(r'Capillary pressure, $Pc [psi]$')
    plt.ylim([-2,100])
    plt.xlim([0,1])
    plt.savefig('capillarypres.png',bbox_inches='tight', dpi = 600)