import matplotlib.pyplot as plt
import numpy as np
from input_file import inputfile

'''
inputfile(fluid,reservoir,numerical,BC,IC)     #uploading all input properties 

def analytical(time,blocks):
    tol = 0.0001
    x = np.arange(blocks)*dx + dx/2
    P1 = 500-((D-reservoir.length*np.sin(dip*np.pi/180))*(rho/144))- 14.7; #based on boundary 500 psi at L P1 at pooint n check
    P0 = 0
    alpha = 6.33E-3*k/(phi*u*ct)
    n = 0
    SUM = 0
    for i in range(10000000):
        new = (-1)**(n+1)*np.cos((n+1/2)*np.pi*x/Length)*np.exp(-(2*n+1)**2/Length**2*alpha*np.pi**2*time/4)/(2*n+1)
        SUM = 4/np.pi*(P1-P0)*new + SUM
        error=max(abs(new));
        if error > tol:
            n = n +1
        else:
            break

    P = SUM + p + P1
    return P;
'''

#fluid, reservoir and simulation parameters   
def post_process(P_plot,xc,time,nmax):

#plotting
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plot = plt.plot(xc,P_plot[:,0],'r-',label=r"$%0.1f$ [days]" % time[0])
    plt.plot(xc,P_plot[:,10],'r--',label=r"$%0.1f$ [days]" % time[10])
    plt.plot(xc,P_plot[:,20],'r-.',label=r"$%0.1f$ [days]" % time[20])
    plt.plot(xc,P_plot[:,30],'b-',label=r"$%0.1f$ [days]" % time[30])
    plt.plot(xc,P_plot[:,40],'b--',label=r"$%0.1f$ [days]" % time[40])
    plt.plot(xc,P_plot[:,50],'b-.',label=r"$%0.1f$ [days]" % time[50])
    plt.plot(xc,P_plot[:,100],'g-',label=r"$%0.1f$ [days]" % time[100])
    plt.plot(xc,P_plot[:,150],'g--',label=r"$%0.1f$ [days]" % time[150])
    plt.plot(xc,P_plot[:,200],'g-.',label=r"$%0.1f$ [days]" % time[200])
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    plt.ylabel(r'$Pressure$ [psi]')
    plt.xlabel(r'$x$ [feet]')
    legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig('Pvsx.png',bbox_inches='tight', dpi = 600)