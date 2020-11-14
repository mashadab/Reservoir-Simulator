'''A short subroutine to solve A*x = B where A is pentadiagonal, square and real matrix and B is a vector
Uses pentapy and numpy libraries | More than 10 times faster than spsolve

Cite: Müller, (2019). pentapy: A Python toolbox for pentadiagonal linear systems. 
Journal of Open Source Software, 4(42), 1759, https://doi.org/10.21105/joss.01759

Author: Mohammad Afzal Shadab
Email : mashadab@utexas.edu
Date  : 10/15/2020
'''

import pentapy as pp
import numpy   as np
from scipy.sparse import issparse
'''
# Function to extract diagonals
def flatten_diagonally(x, diags=None):
    diags = np.array(diags)
    if x.shape[1] > x.shape[0]:
        diags += x.shape[1]-x.shape[0]
    n = max(x.shape)
    ndiags = 2*n-1
    i,j = np.indices(x.shape)
    d = np.array([])
    for ndi in range(ndiags):
        if diags != None:
            if not ndi in diags:
                continue
        d = np.concatenate((d,x[i==j+(n-1)-ndi]))
    return d

# Function to solve Ax = B
def pentapy_solveA(A,B):
    #finding the size of the matrix
    if issparse(A): 
        size,size_dummy = A.shape
    else:
        size = len(A)

    #extracting the diagonals
    second_super = [np.append(flatten_diagonally(A, diags=(size-1+2)),[0.0, 0.0])]  #second subdiagonal
    first_super  = [np.append(flatten_diagonally(A, diags=(size-1+1)),[0.0])]       #first subdiagonal
    main         = [flatten_diagonally(A, diags=(size-1))]                          #main diagonal
    first_sub    = [np.append([0.0],flatten_diagonally(A, diags=(size-1-1)))]       #first superdiagonal
    second_sub   = [np.append([0.0, 0.0],flatten_diagonally(A, diags=(size-1-2)))]  #second superdiagonal
    A_diagflat   = np.concatenate((second_super, first_super,main,first_sub,second_sub), axis=0)
    
    #solving the system of equation
    X = pp.solve(A_diagflat, B, is_flat=True)
    return X;
    
# Function to solve Ax = B
def pentapy_solveB(A,B):
    #finding the size of the matrix
    if issparse(A): 
        size,size_dummy = A.shape
    else:
        size = len(A)

    #extracting the diagonals
    second_super = [np.append(A.diagonal(2),[0.0, 0.0])]  #second subdiagonal
    first_super  = [np.append(A.diagonal(1),[0.0])]       #first subdiagonal
    main         = [A.diagonal()]                          #main diagonal
    first_sub    = [np.append([0.0],A.diagonal(-1))]       #first superdiagonal
    second_sub   = [np.append([0.0, 0.0],A.diagonal(-2))]  #second superdiagonal
    A_diagflat   = np.concatenate((second_super, first_super,main,first_sub,second_sub), axis=0)
    
    #solving the system of equation
    X = pp.solve(A_diagflat, B, is_flat=True)
    return X;
'''

# Function to solve Ax = B
def pentapy_solve(A,B):
    #extracting the diagonals
    second_super = [np.append(A.diagonal(2),[0.0, 0.0])]   #second subdiagonal
    first_super  = [np.append(A.diagonal(1),[0.0])]        #first subdiagonal
    main         = [A.diagonal()]                          #main diagonal
    first_sub    = [np.append([0.0],A.diagonal(-1))]       #first superdiagonal
    second_sub   = [np.append([0.0, 0.0],A.diagonal(-2))]  #second superdiagonal
    A_diagflat   = np.concatenate((second_super, first_super,main,first_sub,second_sub), axis=0)
    
    #solving the system of equation
    X = pp.solve(A_diagflat, B, is_flat=True)
    return X;


