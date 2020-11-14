"""
reservoir simulation assignment 3
1D reservoir simulation Q3.4: 
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 9/24/2020
"""
import numpy as np
import matplotlib.pyplot as plt


def postprocess(P_plot,numerical,well,time):
    nmax = int(numerical.tfinal / numerical.dt)
    Pwf  = np.empty((nmax,len(well.x)))
    qwf  = np.empty((nmax,len(well.x)))

    '''
    for j in range(0,len(well.x)):
        for i in range(0,nmax):
            if well.type[j][0] == 1:  # 1 for rate  
                Pwf[i,j]  = well.constraint[j][0] / well.Jvec[j][0] + P_plot[well.block[j][0],i]
                qwf[i,j]  = well.constraint[j][0]
            elif (well.type[j][0] == 2):# 2 for BHP  
                Pwf[i,j]  = well.constraint[j][0]
                qwf[i,j]  = well.Jvec[j][0] * (P_plot[well.block[j][0],i] - well.constraint[j][0] )
    '''
    for j in range(0,len(well.x)):
        for i in range(0,nmax):
            if well.typetime[i][j] == 1:  # 1 for rate  
                Pwf[i,j]  = well.constrainttime[i][j] / well.Jvec[j][0] + P_plot[well.block[j][0],i]
                qwf[i,j]  = well.constrainttime[i][j]
            elif (well.typetime[i][j] == 2):# 2 for BHP  
                Pwf[i,j]  = well.constrainttime[i][j]
                qwf[i,j]  = well.Jvec[j][0] * (P_plot[well.block[j][0],i] - well.constrainttime[i][j] )    

    fig = plt.figure()
    plot = plt.plot(time[0:len(time)-1],qwf[:,0],'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],qwf[:,1],'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],qwf[:,2]+qwf[:,3]+qwf[:,4],'b',label=f'Well 3')
    plt.plot(time[0:len(time)-1],qwf[:,5]+qwf[:,6]+qwf[:,7],'g--',label=f'Well 4')
    plt.plot(time[0:len(time)-1],qwf[:,8],'c',label=f'Well 5')
    plt.plot(time[0:len(time)-1],qwf[:,9],'m--',label=f'Well 6')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Flow rate $[scf/day]$')
    plt.xlabel(r'time [days]')
    plt.savefig('wellflowratesvstime.png',bbox_inches='tight', dpi = 600)
    
    fig = plt.figure()
    plot = plt.plot(time[0:len(time)-1],Pwf[:,0],'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],Pwf[:,1],'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],Pwf[:,2],'b',label=f'Well 3')
    plt.plot(time[0:len(time)-1],Pwf[:,5],'g--',label=f'Well 4')
    plt.plot(time[0:len(time)-1],Pwf[:,8],'c',label=f'Well 5')
    plt.plot(time[0:len(time)-1],Pwf[:,9],'m--',label=f'Well 6')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'BHP [psi]')
    plt.xlabel(r'time [days]')
    plt.savefig('BHPvstime.png',bbox_inches='tight', dpi = 600)