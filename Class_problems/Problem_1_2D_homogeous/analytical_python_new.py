import numpy as np

def analytical(time,blocks):
    tol = 0.00000001
    rho = 62.4
    Length = 10000
    za   = 1000/(rho/144.0)
    dip = 0.0*np.pi/6.0
    k   = 50.0
    phi = 0.2
    ct  = 1E-6
    mu  = 1.0
    dx = Length / blocks
    x = np.arange(blocks)*dx + dx/2
    P1 = 2000-((za-Length*np.sin(dip))*(rho/144)) - 14.6959502543; #based on boundary 500 psi at L P1 at pooint n check
    P0 = 0
    alpha = 6.33E-3*k/(phi*mu*ct)
    n = 0
    SUM = 0
    for i in range(10000000):
        new = (-1)**(n+1)*np.cos((n+1/2)*np.pi*x/Length)*np.exp(-(2*n+1)**2/Length**2*alpha*np.pi**2*time/4)/(2*n+1)
        SUM = 4/np.pi*(P1-P0)*new + SUM
        error=max(abs(new));
        if error > tol:
            n = n + 1
        else:
            break
    dpg=((za-x*np.sin(dip))*(rho/144.0)) + 14.6959502543
    P = SUM + dpg + P1
    P = np.flip(P)
    return P, x