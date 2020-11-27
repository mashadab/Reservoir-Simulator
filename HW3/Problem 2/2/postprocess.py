
import numpy as np
import  matplotlib.pyplot as plt

Nthree   = np.load('tvsC_N3X3.npz')
Nthirty  = np.load('tvsC_N30X30.npz') 

plt.figure()
plt.plot(Nthree['arr_0'],Nthree['arr_1'],'r-',label='N=3X3')
plt.plot(Nthirty['arr_0'],Nthirty['arr_1'],'k-',label='N=30X30')
plt.legend(loc='best', shadow=False, fontsize='medium')
plt.xlabel(r'$t_D$')
plt.ylabel(r'$C_D$')
plt.savefig('CdvsTd.png',bbox_inches='tight', dpi = 600)