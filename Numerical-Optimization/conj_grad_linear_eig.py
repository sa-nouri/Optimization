# -*- coding: utf-8 -*-

import random
import numpy as np
import matplotlib.pyplot as plt

# Linear Equation Ax=b

rk = 1e-6

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

n = 50
x = np.zeros((n, 1))
b = np.ones(shape=(n, 1))
eigenvalues1 = np.random.uniform(low=1.0, high=8.0, size=n)

centers = [2, 5, 10]
sigma = [0.35, 0.65, 1.0]
sample1 = np.random.normal(centers[0], sigma[0], int(n/3))
sample2 = np.random.normal(centers[1], sigma[1], int(n/3))
sample3 = np.random.normal(centers[2], sigma[2], n - 2*int(n/3))

data = []
data = data + list(sample1) + list(sample2) + list(sample3)
eigenvalues2 = random.sample(data, n)

A1 = np.diag(eigenvalues1)
A2 = np.diag(eigenvalues2)
print(np.linalg.cond(A1))
print(np.linalg.cond(A2))

xs1, curve_x1 = conjgrad(A1, b, x)
xs2, curve_x2 = conjgrad(A2, b, x)

distance1 = curve_x1 - xs1
distance2 = curve_x2 - xs2

result1 = np.linalg.norm(distance1, ord=1, axis=1, keepdims=False) ** 2
result2 = np.linalg.norm(distance2, ord=1, axis=1, keepdims=False) ** 2

plt.figure(figsize=(8, 8))
plt.plot(np.log(result1), '-b')
plt.plot(np.log(result2), '-r')
plt.grid()
plt.legend(['Uniform Eigenvalues', 'Clustered Eigenvalues'])
plt.xlabel('Iterations')
plt.ylabel('logarithm')
