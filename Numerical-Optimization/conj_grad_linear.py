# -*- coding: utf-8 -*-
import numpy as np

def conjgrad(A, b, x):
    """
    A function to solve [A]{x} = {b} linear equation system with the 
    conjugate gradient method.
    ========== Parameters ==========
    A : matrix 
        A real symmetric positive definite matrix.
    b : vector
        The right hand side (RHS) vector of the system.
    x : vector
        The starting guess for the solution.
    """  
    r = b - np.dot(A, x)
    p = r
    rsold = np.dot(r.T, r)
    curve_x = []
    while np.sqrt(rsold) > rk:
        Ap = np.dot(A, p)
        alpha = rsold / np.dot(p.T, Ap)
        x = x + (alpha * p)
        r = r - (alpha * Ap)
        rsnew = np.dot(r.T, r)
        p = r + (rsnew/rsold)*p
        rsold = rsnew
        curve_x.append(list(x))
    return x, np.array(curve_x)

# Linear Equation Ax=b

n = 20 # 5, 8, 12, 20
rk = 1e-6

x = np.zeros((n, 1))
b = np.ones(shape=(n, 1))
A = np.zeros((n, n))

for i in range(1, n+1):
    for j in range(1, n+1):
        A[i-1, j-1] = 1/(i + j -1)

# xs, array = conjgrad(A, b, x)

array.shape
np.linalg.cond(A)