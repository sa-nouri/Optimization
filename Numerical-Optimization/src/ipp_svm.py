# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import collections
import sys
import scipy
import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.model_selection import train_test_split

def get_data():
    np.random.seed(27)
    X, y = datasets.fetch_openml('mnist_784', version=1, return_X_y=True)
    sub_indices = (y == '1') | (y == '0')
    X = X[sub_indices]
    y = y[sub_indices]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
    return X_train, X_test, y_train, y_test

whole_X_train, X_test, whole_y_train, y_test = get_data()
# split further into train and val 
X_train, X_val, y_train, y_val = train_test_split(whole_X_train, whole_y_train, 
                                                 test_size=0.15, random_state=42)

y_train_floats = [float(i) for i in y_train]
y_val_floats = [float(i) for i in y_val]
y_test_floats = [float(i) for i in y_test]
whole_y_train_floats = [float(i) for i in whole_y_train]

label_map = {}
unique_labels = set(y_train_floats)
assert(len(unique_labels) == 2)
flag = True       
for label in unique_labels:
    if flag:
        label_map[label] = -1
        flag = False
    else:
        label_map[label] = 1
print(label_map)

y_train_floats = np.array([label_map[label] for label in y_train_floats])
y_val_floats = np.array([label_map[label] for label in y_val_floats])
y_test_floats = np.array([label_map[label] for label in y_test_floats])
whole_y_train_floats = np.array([label_map[label] for label in whole_y_train_floats])

class SVMObjective:
    """
    Implements functions and derivatives for SVM Objective. Note that the
    function is the barrier modified function, not the original dual SVM
    objective
    """
    def __init__(self, X_train, y_train, kernel_func, C):
       self.K = kernel_func(X_train, X_train)
       self.Y = np.diag(y_train)
       self.C = C

    def f(self, a, t):
        """
        Evaluate function f
        """
        assert a.shape == (self.Y.shape[0], 1)
        first =  .5 * a.T @ self.Y @ self.K @ self.Y @ a
        second = np.sum(a)
        third = np.sum(np.log(a))
        fourth = np.sum(np.log(self.C - a))
        return t * (first - second) - third - fourth

    def grad(self, a, t):
        """
        Returns gradient of SVM barrier function
        """
        first = t * (self.Y @ self.K @ self.Y @ a  - np.ones_like(a))
        second = np.array(1.0/a)
        third = np.array(1.0/(self.C - a))

        return first - second + third

    def hessian(self, a, t):
        """
        Returns hessian of SVM barrier function
        """
        first = t * self.Y @ self.K @ self.Y
        second = np.diag(1/a**2)
        third = np.diag(1/(self.C - a)**2)
        return first + second + third

class Newton:
    def __init__(self, curr_a, objective, C, c_1=0.01, rho=0.5):
        self.curr_a = curr_a
        self.objective = objective
        self.c_1 =c_1
        self.rho = rho
        self.A = np.vstack((np.eye(curr_a.shape[0]),
            -1*np.eye(curr_a.shape[0])))
        self.C = C

    def newton_step(self, t, newton_decrement=False):
        hess = self.objective.hessian(self.curr_a, t)
        grad = self.objective.grad(self.curr_a, t)
        KKT_mat = np.block([[hess, self.A.T], [self.A,
            np.zeros((self.A.shape[0], self.A.T.shape[1]))]])
        rhs = -1 * np.vstack((grad, np.zeros((self.A.shape[0], 1))))

        try:
            solved = scipy.linalg.solve(KKT_mat, rhs)
        except np.linalg.LinAlgError:
            solved = scipy.linalg.pinv(KKT_mat) @ rhs
        delta_a = solved[:hess.shape[0]] # grab just delta x part

        # we calcualte the newton decrement, -grad^T * delta_a, as a stopping
        # criterion
        newton_decrement = -1 * np.dot(grad.T, delta_a)
        return delta_a, newton_decrement

    def line_search(self, delta_x, t):
        """
        Perform line search on Newton direction from delta_x
        Backtracks until sufficient decrease (Armijo) cond is satisfied and the inequality
        constraints are also satisfied
        Returns step size
        """

        ls_coeff = 1
        # Inequality constraint condition
        while (min(self.curr_a + ls_coeff * delta_x) < 0 or max(self.curr_a +
            ls_coeff * delta_x) > self.C):
            ls_coeff = self.rho * ls_coeff

        # Sufficient decrease conditions
        while True:
            lower_bound = self.objective.f(self.curr_a, t) + self.c_1 * ls_coeff * np.dot(self.objective.grad(self.curr_a, t).T, delta_x)
            if self.objective.f(self.curr_a + ls_coeff * delta_x, t) <= lower_bound:
                break
            ls_coeff = self.rho * ls_coeff

        return ls_coeff

class Barrier:
    """
    Implements Barrier method (i.e. outer loop logic)
    """
    def __init__(self, t_0=1, mu=3, tol=1e-10, newton_tol = 1e-3, max_newton_iter=10, max_iter=200):
        self.t_0 = t_0
        self.mu = mu # expansion constant
        self.tol = tol
        self.newton_tol = newton_tol
        self.max_newton_iter = max_newton_iter
        self.max_iter = max_iter

    def run(self, X_train, y_train, kernel_func, C, m, store_iterates=False):
        """
        Runs Barrier method
        Parameters:
        X_train: data as rows
        y_train: array of floats representing labels, +1 and -1
        kernel_func: callable kernel function
        C: SVM C parameter
        m: number of inequality constraints
        """
        a = self._feasible_starting_point(y_train, C)
        svm_objective = SVMObjective(X_train, y_train, kernel_func, C)
        t = self.t_0
        curr_iter = 0
        if store_iterates:
            iterates = [a]
        while True:
            """
            Includes newton decrement for checking inner loop Newton
            convergence, *although* it is enough to solve the inner loop
            imperfectly because we only care about the centering path, not
            the solutions on the path itself
            """
            curr_newton_iter = 0
            while True:
                newton = Newton(a, svm_objective, C)
                unscaled_da, newton_decrement = newton.newton_step(t)
                step_length = newton.line_search(unscaled_da, t)
                a_next = a + step_length * unscaled_da
                a_old = a
                a = a_next
                """
                Check stopping criteria for inner loop Newton.
                Empirically sometimes the step is incredibly tiny but the
                Newton decrement is larger than the default tol of 1e-3 -- the
                np.isclose condition and the max Newton iterations
                """
                if newton_decrement < self.newton_tol or np.all(np.isclose(a, a_old)) or curr_newton_iter >= self.max_newton_iter:
                    break
                curr_newton_iter +=1
            # increase t
            t_new = t * self.mu
            print("Finished iter: {}".format(curr_iter))
            curr_iter+=1
            if store_iterates:
                iterates.append(a)

            # check stopping criterion for outer loop barrier
            if m/t < self.tol:
                break
            t = t_new
        if store_iterates: return a, iterates
        return a

    def _feasible_starting_point(self, y_labels, C):
        """
        Find a strictly feasible starting point for the barrier method
        """
        total_num = len(y_labels)
        pos = np.sum(y_labels == 1)
        pos_frac = pos / (total_num - pos)
        a0 = np.zeros((y_labels.shape[0], 1))
        for i in range(a0.shape[0]):
            if y_labels[i] == 1:
                a0[i] = C * (1 - pos/total_num)
            else:
                a0[i] = C * pos_frac * (1 - pos/total_num)
        return a0

class SVMClassifier:
    """
    Takes solved Lagrange multipliers barrier method and creates SVM classifier
    """

    def __init__(self, a, X_train, y_train, kernel_func):
        self.alphas = a
        if len(y_train.shape) != 2 or y_train.shape[1] != 1: y_train = y_train[:,
                None]
        alpha_y = self.alphas * y_train
        if len(alpha_y.shape) != 2 or alpha_y.shape[1] != 1: alpha_y = alpha_y[:, None]
        self.alpha_y = alpha_y
        # calculate b
        K_train = kernel_func(X_train, X_train)
        support_vectors = np.where(self.alphas != 0)[0]
        differences = []
        for ind in support_vectors:
            kernels = K_train[ind][:, None] # kernel value of support vector w train pts
            pred = kernels.T @ self.alpha_y
            diff = np.abs(y_train[ind] - pred)
            differences.append(diff)
        # set b to the median of this differences array
        self.b = np.median(np.array(differences))
        self.kernel_func = kernel_func
        self.X_train = X_train

    def accuracy(self, X_test, y_test, b = None):
        kernel_mat = self.kernel_func(X_test, self.X_train)
        if b is None:
            b = self.b
        pred = np.sign(kernel_mat @ self.alpha_y + b)
        pred = np.squeeze(pred)
        num_correct = sum(pred == y_test)
        print("predictions: {}".format(pred))
        print("true labels: {}".format(y_test))
        print("num correct {}, accuracy {}".format(num_correct,
            num_correct/len(y_test)))

        return num_correct/len(y_test)

def kernel_rbf(gamma):
    def _kernel_rbf(x1, x2):
        if len(x1.shape) == 1:
            x1 = x1[None,:]
        if len(x2.shape) == 1:
            x2 = x2[None, :]
        x_input_norms = np.diag(x1 @ x1.T)[:, None]
        x_out_norms = np.diag(x2 @ x2.T)[None, :]

        numerator = x_input_norms + x_out_norms
        numerator -= 2*(x1 @ x2.T)
        ret = np.exp(-gamma * numerator)

        return ret
    return _kernel_rbf

def kernel_poly(d=3):
    def _kernel_poly(x1, x2):

        if len(x1.shape) == 1:
            x1 = x1[None,:]
        if len(x2.shape) == 1:
            x2 = x2[None, :]
        K = (x1 @ x2.T)**d
        return K
    return _kernel_poly

def kernel_lin():
    def _kernel_linear(x1, x2):
        if len(x1.shape) == 1:
            x1 = x1[None, :]
        if len(x2.shape) == 1:
            x2 = x2[None, :]
        K = x1 @ x2.T
        return K
    return _kernel_linear

barrier = Barrier()
# Cross Validation

best_acc = float('-inf')
best_g = None
best_C = None 

best_accs = []
best_gs = []
best_Cs = []

m = 2 * y_train_floats.shape[0] # number of constraints, two per Lagrange multiplier 
for C in np.logspace(-1, 3, num=5):
    for g in np.logspace(-10, -6, num = 15):
        kernel = kernel_rbf(g)
        lagrange_multipliers = barrier.run(np.array(X_train), y_train_floats, kernel, C, m)
        svm = SVMClassifier(lagrange_multipliers, X_train, y_train_floats, kernel)
        
        # test on acc 
        acc = svm.accuracy(X_val, y_val_floats)
        
        if acc > best_acc:
            best_acc = acc
            best_g = g
            best_C = C 
            
        if acc >= 0.9:
            best_accs.append(acc)
            best_gs.append(g)
            best_Cs.append(C)