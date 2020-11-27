#courant number evaluation   
def courant(i,j,P,numerical,fluid,reservoir):
    
    #calculate velocity using Darcy's law
    u = (reservoir.perm[i,0]/fluid.mu[i,0])*6.33e-3*(P[i]-P[j])/numerical.dx[i,0]   #Darcy velocity ft/day
    
    #calculating interblock Courant number
    cr= abs(u)*numerical.dt / (reservoir.phi[i,0]*numerical.dx[i,0])
    #cr= 0.5
    
    return cr;