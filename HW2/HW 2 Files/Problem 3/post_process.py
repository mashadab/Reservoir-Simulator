import matplotlib.pyplot as plt
from analytical_python_new import analytical

def post_process(P_plot,xc,time,nmax):
    P_plot_analytical_10,x  = analytical(time[10],nmax)
    P_plot_analytical_100,x = analytical(time[100],nmax)
    P_plot_analytical_200,x = analytical(time[200],nmax)
    
    #plotting
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plot = plt.plot(xc,P_plot[:,0],'k-',label=r"C: $%0.1f$ [days]" % time[0])
    plt.plot(xc,P_plot[:,1],'c-',label=r"C: $%0.1f$ [days]" % time[1]) 
    plt.plot(xc,P_plot[:,2],'m-',label=r"C: $%0.1f$ [days]" % time[2]) 
    plt.plot(xc,P_plot[:,5],'y-',label=r"C: $%0.1f$ [days]" % time[5]) 
    plt.plot(xc,P_plot[:,10],'r-',label=r"C: $%0.1f$ [days]" % time[10])    
    #plt.plot(x,P_plot_analytical_10,'k--',label=r"A: $%0.1f$ [days]" % time[10])
    plt.plot(xc,P_plot[:,100],'b-',label=r"C: $%0.1f$ [days]" % time[100])
    #plt.plot(x,P_plot_analytical_100,'k--',label=r"A: $%0.1f$ [days]" % time[100])   
    plt.plot(xc,P_plot[:,200],'g-',label=r"C: $%0.1f$ [days]" % time[200])
    #plt.plot(x,P_plot_analytical_200,'k--',label=r"A: $%0.1f$ [days]" % time[200])     
    #manager = plt.get_current_fig_manager()
    #manager.window.showMaximized()
    plt.ylabel(r'$Pressure$ [psi]')
    plt.xlabel(r'$x$ [feet]')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig('Pvsx.png',bbox_inches='tight', dpi = 600)