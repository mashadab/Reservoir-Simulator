'''A short subroutine to solve A*x = B where A is pentadiagonal, square and real matrix and B is a vector
Uses pentapy and numpy libraries | More than 100 times faster than spsolve

Cite: Müller, (2019). pentapy: A Python toolbox for pentadiagonal linear systems. 
Journal of Open Source Software, 4(42), 1759, https://doi.org/10.21105/joss.01759

Author: Mohammad Afzal Shadab
Email : mashadab@utexas.edu
Date  : 10/15/2020
'''

import numpy as np
import pentapy as pp
from scipy.sparse import csr_matrix
from pentapy_solve import pentapy_solve
from scipy.sparse.linalg import spsolve
import time as timer

size = 10000
# create a flattened pentadiagonal matrix
M_flat = (np.random.random((5, size)) - 0.5) * 1e-5
M = pp.create_full(M_flat, col_wise=False)
V = np.random.random(size) * 1e5

Msparse = csr_matrix(M)

penta_start = timer.clock()
X_new = pentapy_solve(Msparse,V)
penta_end = timer.clock()
print('Xnew',X_new)

#sparse solve no magic
spsolve_start = timer.clock()
X_spsolve = spsolve(Msparse,V)
spsolve_end = timer.clock()
print('Xspsolve',X_spsolve)

# calculate the error
print(np.max(np.abs(np.dot(M, X_new) - V)))
print(np.max(np.abs(np.dot(M, X_spsolve) - V)))
print('Time, Pentapy:',spsolve_end-spsolve_start,'(s), spsolve:',penta_end-penta_start,'(s)')