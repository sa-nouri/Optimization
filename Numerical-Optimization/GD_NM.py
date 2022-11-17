#@title Import libraries

import numpy as np
import numpy.linalg as npla
import matplotlib.pyplot as plt

from numpy import linalg as LA

#@title Define function, gradient, and hessian

def f(x):
    return (1 - x[0])**2+ 100*((x[1]-x[0]**2)**2)

def f_dx(x):
    df1 = -2*(1 - x[0]) - (400*x[0])*(x[1] - (x[0]**2))
    df2 = 200*(x[1] - (x[0]**2))
    return np.array([df1, df2])

def f_hessian(x):
    d2f_dx2 = 2 - 400*x[1] + 1200 * (x[0]**2)
    d2f_dyx = -400*x[0]
    d2f_dy2 = 200
    return(np.matrix([[d2f_dx2, d2f_dyx], [d2f_dyx, d2f_dy2]]))

#@title LINE SEARCH STEP SIZE

def backtrack2(x0, f, fdx, t = 1, alpha=0.2, beta=0.8):
    while f(x0 - t*fdx(x0)) > f(x0) + alpha * t * np.dot(fdx(x0).T, -1*fdx(x0)):
        t *= beta
    return t

#@title Steepest Descent

def grad(point, max_iter, f, fdx):
    iter = 1
    points = []
    values = []
    step_sizes = []
    while (np.linalg.norm(fdx(point)) > 0.000001):

        t = backtrack2(point, f, fdx)
        point = point - np.dot(t, fdx(point))
        iter += 1
        if iter > max_iter:
            break
        points.append(point)
        values.append(f(point))
        step_sizes.append(t)

    return points, values, step_sizes

#@title Newton's Method

def lambda_sq(fdx, Hessian, point):
    val_hessian = Hessian(point)
    if val_hessian.ndim < 2:
        lambda_sq = fdx(point) * val_hessian * fdx(point)
    else:
        lambda_sq = np.dot(np.dot(fdx(point), npla.pinv(val_hessian)), fdx(point).T)
    return (np.asscalar(lambda_sq))

def delta_x(fdx, Hessian, point):
    val_hessian = Hessian(point)
    if val_hessian.ndim < 2:
        delta_x = val_hessian * fdx(point)
    else:
        delta_x = np.dot(-npla.pinv(Hessian(point)), fdx(point).T)
    return (delta_x)

def newtons_method(x, f, fdx, Hessian, eps=0.00001, max_iters=20000):
    iters = 1
    points = []
    values = []
    step_sizes = []
    lmb_sq = lambda_sq(fdx, Hessian, x)
    while ((lmb_sq/2.0) > eps):
        dlt_x = delta_x(fdx, Hessian, x)
        t = backtrack2(x, f, fdx)
        x = np.array((x + np.dot(t , dlt_x)))[0]
        lmb_sq = lambda_sq(fdx, Hessian, x)
        
        iters += 1
        if(iters > max_iters):
            break
            
        points.append(x)
        values.append(f(x))
        step_sizes.append(t)
    return points, values, step_sizes

#@title Simulation

point1 = [6.0, 6.0]
point2 = [2.0, 0.0]
point3 = [5.0, -6.0]
point4 = [-1.0, 4.0]
points, values, step_sizes = grad(point1, max_iter=20000, f=f, fdx=f_dx)
points1, values1, step_sizes1 = grad(point2, max_iter=20000, f=f, fdx=f_dx)
points2, values2, step_sizes2 = grad(point3, max_iter=20000, f=f, fdx=f_dx)
points3, values3, step_sizes3 = grad(point4, max_iter=20000, f=f, fdx=f_dx)
# points, values, step_sizes = newtons_method(point, f, f_dx, f_hessian,
                                            # eps=0.00001, max_iters=20000)

#@title Plot function values

plt.figure()
plt.plot(values, '-b')
plt.plot(values1, '-r')
plt.plot(values2, '-g')
plt.plot(values3, '-y')
plt.ylabel('Function Values')
plt.xlabel('Iterations')
plt.grid()
# plt.title("Newton Method")
plt.title("Steepest Descent Method")

#@title plot distance to distance to optimal point

points1 = points1 - points1[-1]
dist1 = LA.norm(points1, axis=1)

points2 = points2 - points2[-1]
dist2 = LA.norm(points2, axis=1)

points3 = points3 - points3[-1]
dist3 = LA.norm(points3, axis=1)

points = points - points[-1]
dist = LA.norm(points, axis=1)


plt.figure()
plt.plot(dist, 'b')
plt.plot(dist1, 'r')
plt.plot(dist2, 'g')
plt.plot(dist3, 'y')
plt.ylabel('Distance to optimal point')
plt.xlabel('Iterations')
plt.grid()
# plt.title("Newton Method")
plt.title("Steepest Descent Method")

#@title plot distance to step sizes

plt.figure()
plt.plot(step_sizes, 'b')
plt.plot(step_sizes1, 'r')
plt.plot(step_sizes2, 'g')
plt.plot(step_sizes3, 'y')
plt.ylabel('Step size')
plt.xlabel('Iterations')
plt.grid()
# plt.title("Newton Method")
plt.title("Steepest Descent Method")