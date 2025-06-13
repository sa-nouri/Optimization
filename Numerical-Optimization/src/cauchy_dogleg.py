# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as ln
import scipy as sp
import matplotlib.pyplot as plt

from math import sqrt
from numpy import linalg as la


def is_psd(a):
    if a.shape == (1,):
        return a >= 0
    else:
        return np.all(la.eigvals(a) >= 0)


def repair_psd(a, eps=1e-3):
    if is_psd(a):
        return a
    if a.shape == (1,):
        return np.atleast_1d(eps)
    else:
        min_eigval = min(la.eigvals(a))
        if min_eigval >= 0:
            return a
        else:
            x = a - (1+eps)*min_eigval*np.identity(a.shape[0])
            return x


def inverse_matrix(a):
    if a.shape == (1,):
        return np.atleast_1d(1. / a)
    else:
        return la.inv(a)

# Objective function
def f(x):
    return np.array((x[1] - 0.129*(x[0]**2) + 1.6*x[0]-6)**2 + 6.07*np.cos(x[0]) + 10)


# Gradient
def jac(x):
    return np.array([2*(-0.129*x[0]*2+1.6) * (x[1] - 0.129*x[0]**2 + 1.6*x[0]-6) - np.sin(x[0]),
                    2*(x[1] - 0.129*x[0]**2 + 1.6*x[0]-6)])


# Hessian
def hess(x):
    a = -2*0.129*2*(x[1] - 0.129*x[0]**2 + 1.6*x[0]-6) + 2*(-0.129*x[0]*2+1.6) * (2*0.129*x[0] + 1.6) - np.cos(x[0])
    b = 2*(-0.129*x[0]*2+1.6)
    c = -2*0.129*2*x[0] + 3.2
    d = 2
    hess = [[a, b], [c, d]]
    return np.array(hess)

def dogleg_method(Hk, gk, Bk, trust_radius):

    # Compute the Newton point.
    # This is the optimum for the quadratic model function.
    # If it is inside the trust radius then return this point.
    pB = -np.dot(Hk, gk)
    norm_pB = sqrt(np.dot(pB, pB))

    # Test if the full step is within the trust region.
    if norm_pB <= trust_radius:
        return pB

    # Compute the Cauchy point.
    # This is the predicted optimum along the direction of steepest descent.
    pU = - (np.dot(gk, gk) / np.dot(gk, np.dot(Bk, gk))) * gk
    dot_pU = np.dot(pU, pU)
    norm_pU = sqrt(dot_pU)

    # If the Cauchy point is outside the trust region,
    # then return the point where the path intersects the boundary.
    if norm_pU >= trust_radius:
        return trust_radius * pU / norm_pU

    pB_pU = pB - pU
    dot_pB_pU = np.dot(pB_pU, pB_pU)
    dot_pU_pB_pU = np.dot(pU, pB_pU)
    fact = dot_pU_pB_pU ** 2 - dot_pB_pU * (dot_pU - trust_radius ** 2)
    tau = (-dot_pU_pB_pU + sqrt(fact)) / dot_pB_pU

    # Decide on which part of the trajectory to take.
    return pU + tau * pB_pU

def cauchy_point(x, b_matrix, delta):
    g_matrix = f_provider.grad(*x)
    gT_b_g = np.dot(np.dot(g_matrix, b_matrix), g_matrix)
    gT_g = np.dot(g_matrix, g_matrix)
    g_norm = norm_p(g_matrix)

    if gT_b_g > 0 and abs(gT_g / gT_b_g) * g_norm < delta:
        alpha = gT_g / gT_b_g
    else:
        alpha = delta / g_norm

    return -alpha * g_matrix


def trust_region_dogleg(func, jac, hess, x0, initial_trust_radius=1.0,
                        max_trust_radius=100.0, eta=0.15, gtol=1e-4,
                        maxiter=100):
    xk = x0
    trust_radius = initial_trust_radius
    k = 0
    while True:

        gk = jac(xk)
        Bk = hess(xk)
        Hk = np.linalg.inv(Bk)

        pk = dogleg_method(Hk, gk, Bk, trust_radius)

        # Actual reduction.
        act_red = func(xk) - func(xk + pk)

        # Predicted reduction.
        pred_red = -(np.dot(gk, pk) + 0.5 * np.dot(pk, np.dot(Bk, pk)))

        # Rho.
        rhok = act_red / pred_red
        if pred_red == 0.0:
            rhok = 1e99
        else:
            rhok = act_red / pred_red

        # Calculate the Euclidean norm of pk.
        norm_pk = sqrt(np.dot(pk, pk))

        # Rho is close to zero or negative, therefore the trust region is shrunk.
        if rhok < 0.25:
            trust_radius = 0.25 * norm_pk
        else:
            # Rho is close to one and pk has reached the boundary of the trust region, therefore the trust region is expanded.
            if rhok > 0.75 and norm_pk == trust_radius:
                trust_radius = min(2.0 * trust_radius, max_trust_radius)
            else:
                trust_radius = trust_radius

        # Choose the position for the next iteration.
        if rhok > eta:
            xk = xk + pk
        else:
            xk = xk

        # Check if the gradient is small enough to stop
        if ln.norm(gk) < gtol:
            break

        # Check if we have looked at enough iterations
        if k >= maxiter:
            break
        k = k + 1
    return xk


result = trust_region_dogleg(f, jac, hess, [6, 14])
print("Result of trust region dogleg method: {}".format(result))
print("Value of function at a point: {}".format(f(result)))

def cauchy_point_step_finder(gx, b, delta):
    gt_b_g = np.matmul(np.atleast_1d(np.matmul(gx.T, b)), gx)
    g_norm = la.norm(gx)
    if gt_b_g <= 0:
        taw = 1
    else:
        taw = min(g_norm**3./(delta*gt_b_g), 1.)
    cp = -1. * (taw*delta/g_norm) * gx
    mul = np.floor(delta / la.norm(cp).astype('f')) * (1-1e-3)
    if mul == 0:
        return cp * (1-1e-3)
    else:
        return mul * cp


def solve_taw_for_dogleg(pu, pb, delta):
    a = np.dot(pu-pb, pu-pb)
    b = -2 * (2*np.dot(pu, pu) + np.dot(pb, pb) - 3*np.dot(pu, pb))
    c = np.dot(2*pu-pb, 2*pu-pb) - delta**2
    d = np.sqrt(b**2 - 4*a*c)
    t1 = (-b + d) / (2*a)
    t2 = (-b - d) / (2*a)
    if 0 <= t1 <= 2:
        if 0 <= t2 <= 2:
            return min(t1, t2)
        return t1
    elif 0 <= t2 <= 2:
        return t2
    else:
        raise ArithmeticError('Taw is not in [0,2]: %d %d', t1, t2)


def dogleg_step_finder(gx, b, delta):
    pb = -np.atleast_1d(np.matmul(inverse_matrix(b), gx))
    if la.norm(pb) <= delta:
        return pb
    pu = - (np.matmul(gx.T, gx).astype('f') / (np.matmul(np.atleast_1d(np.matmul(gx.T, b)), gx))) * gx
    if la.norm(pu) >= delta:
        return (delta / la.norm(pu).astype('f')) * (1-1e-3) * pu
    taw = solve_taw_for_dogleg(pu, pb, delta)
    if taw <= 1:
        return taw * pu
    else:
        return pu + (taw - 1)*(pb - pu)

def model(f, g, b, x, p, delta):
    if la.norm(p) > delta:
        raise ArithmeticError('P must not be bigger than delta, p and delta:', p, la.norm(p), delta)
    return f(x) + np.matmul(g(x).T, p) + .5*np.matmul(np.atleast_1d(np.matmul(p.T, b)), p)


def trust_region(f, g, hf, x0, delta_0, max_delta, etha, step_finder, repair_hessian=True, eps=1e-5):

    x = x0
    delta = delta_0
    iterations = 0

    while True:

        iterations += 1
        b = hf(x)
        if repair_hessian:
            b = repair_psd(b)

        p = step_finder(g(x), b, delta)
        rho = (f(x) - f(x+p)).astype('f') / (model(f, g, b, x, p, delta) - model(f, g, b, x+p, p, delta))

        if rho < .25:
            delta = .25 * delta
        elif rho >= .75 and np.isclose(la.norm(p), delta, 1e-4):
            delta = min(2*delta, max_delta)

        if rho > etha:
            x = x + p
        elif np.allclose(p, np.zeros(p.shape), eps):
            result = x + p
            break

    return result, f(result), iterations

def get_trust_region_information(f, g, hf, x0, delta_0, max_delta, etha, step_finder, repair_hessian=True, eps=1e-5):
    
    x = x0
    delta = delta_0
    iterations = 0
    rhos = []
    original_bs = []

    while True:

        iterations += 1
        b = hf(x)
        if is_psd(b):
            original_bs.append(1)
        else:
            original_bs.append(-1)
        if repair_hessian:
            b = repair_psd(b)

        p = step_finder(g(x), b, delta)
        rho = (f(x) - f(x+p)).astype('f') / (model(f, g, b, x, p, delta) - model(f, g, b, x+p, p, delta))
        rhos.append(rho)

        if rho < .25:
            delta = .25 * delta
        elif rho >= .75 and np.isclose(la.norm(p), delta, 1e-4):
            delta = min(2*delta, max_delta)

        if rho > etha:
            x = x + p
        elif np.allclose(p, np.zeros(p.shape), eps):
            break
            
    return rhos, original_bs

xls = [6, 14]
x1s = [np.array(xls), np.array([2, 2,]), np.array([4, 9,])]
# cauchy_point_step_finder, dogleg_step_finder
trust_region_info = [get_trust_region_information(f, jac, hess, x1s[i], 2, 5, .15, step_finder=cauchy_point_step_finder) for i in range(len(x1s))]

plt.rcParams['figure.figsize'] = [17, 17]
for i in range(len(x1s)):
    
    rh = trust_region_info[i][0]
    bs = trust_region_info[i][1]
    x = [j for j in range(len(rh))]
    
    plt.subplot(len(x1s), 2, (i+1)*2-1)
    plt.scatter(x, rh, color='blue', label='Rho')
    plt.xlabel('Iteration (k) for x0: %s' % x1s[i])
    plt.ylabel('Rho for iteration k')
    
    plt.subplot(len(x1s), 2, (i+1)*2)
    plt.scatter(x, bs, color='red', label='Original Bk is psd')
    plt.xlabel('Iteration (k)')
    plt.ylabel('Original Bk is psd(1) or not(-1) for x0: %s' % x1s[i])