"""
reservoir simulation project 1 (2020): Problem 2
2D reservoir simulation: Post processing
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 10/31/2020
"""

import numpy as np
import matplotlib.pyplot as plt

def postprocess(P_plot,numerical,well,time):
    nmax = int(numerical.tfinal / numerical.dt)
    Pwf  = np.zeros((nmax+1,len(well.x)))
    qwf  = np.zeros((nmax+1,len(well.x)))
    
    for j in range(0,len(well.x)):
        for i in range(0,nmax):
            if well.typetime[i][j] == 1:  # 1 for rate  
                if well.constrainttime[i][j] > 0.0:  #injector well
                    well.fw[j][i] = 1.0
                Pwf[i,j]  = P_plot[well.block[j][0],i] + (1 - well.fw[j][i]) * well.constrainttime[i][j] / (well.Jovectime[j][i])
                qwf[i,j]  = well.constrainttime[i][j]
            elif (well.typetime[i][j] == 2):# 2 for BHP  
                Pwf[i,j]  = well.constrainttime[i][j]
                qwf[i,j]  = well.Qwf[j,i]

    #print('BHP at 100 days', Pwf[100,:])
    #print('BHP at 500 days', Pwf[500,:])

    #print('Flow rate at 100 days', qwf[100,:])
    #print('Flow rate at 500 days', qwf[500,:])
    
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plot = plt.plot(time[0:len(time)-1],qwf[0:len(time)-1,0]/5.615,'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],(qwf[0:len(time)-1,1])/5.615,'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],(qwf[0:len(time)-1,2]+qwf[0:len(time)-1,3]+qwf[0:len(time)-1,4])/5.615,'b',label=f'Well 3')
    plt.plot(time[0:len(time)-1],(qwf[0:len(time)-1,5]+qwf[0:len(time)-1,6]+qwf[0:len(time)-1,7])/5.615,'g--',label=f'Well 4')
    plt.plot(time[0:len(time)-1],(qwf[0:len(time)-1,8])/5.615,'c',label=f'Well 5')
    plt.plot(time[0:len(time)-1],(qwf[0:len(time)-1,9])/5.615,'m--',label=f'Well 6')
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Flow rate $[bbl/day]$')
    plt.xlabel(r'time [days]')
    plt.savefig('wellflowratesvstime.png',bbox_inches='tight')
    
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plot = plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,0],'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,1],'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,2],'b',label=f'Well 3')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,5],'g--',label=f'Well 4')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,8],'c',label=f'Well 5')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,9],'m--',label=f'Well 6')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'BHP [psi]')
    plt.xlabel(r'time [days]')
    plt.savefig('BHPvstime.png',bbox_inches='tight')
    
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plot = plt.plot(time[0:len(time)-1],well.fw[0,0:len(time)-1],'k-',label=f'Well 1')
    plt.plot(time[0:len(time)-1],well.fw[1,0:len(time)-1],'r-',label=f'Well 2')
    plt.plot(time[0:len(time)-1],well.fw[2,0:len(time)-1],'b-',label=f'Well 3')
    plt.plot(time[0:len(time)-1],well.fw[5,0:len(time)-1],'g--',label=f'Well 4')
    plt.plot(time[0:len(time)-1],well.fw[8,0:len(time)-1],'c-',label=f'Well 5')
    plt.plot(time[0:len(time)-1],well.fw[9,0:len(time)-1],'m--',label=f'Well 6')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Water cut $f_w$')
    plt.xlabel(r'time [days]')
    plt.savefig('watercutvstime.png',bbox_inches='tight')
    
    return Pwf;