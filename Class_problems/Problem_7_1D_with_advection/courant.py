#courant number evaluation   
def courant(i,j,P,numerical,fluid,reservoir):
    
    #u = (reservoir.perm/fluid.visc)*6.33e-3*(P[i]-P[j])/numerical.dx   #Darcy velocity ft/day
    
    #cr= abs(u)*numerical.dt / (reservoir.phi*numerical.dx)
    
    cr= 0.5
    
    return cr;