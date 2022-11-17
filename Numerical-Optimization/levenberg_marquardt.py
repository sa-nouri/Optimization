# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import collections
import sys
import matplotlib.pyplot as plt

df_train = pd.read_csv('/content/sample_data/mnist_train_small.csv')
df_test = pd.read_csv('/content/sample_data/mnist_test.csv')

train_labels = np.array(df_train['6'])
df_train.drop(columns=['6'], inplace=True)
train_data = np.array(df_train)

test_labels = np.array(df_test['7'])
df_test.drop(columns=['7'], inplace=True)
test_data = np.array(df_test)

# y = f(X, theta) + eps

class Prior:
    
    GAUSSIAN  = 'GAUSSIAN'
    LOGNORMAL = 'LOGNORMAL'
    
    def __init__(self, density, mu = 0., bta = 1.):
        
        assert density in [Prior.GAUSSIAN, Prior.LOGNORMAL], \
               "Invalid density {}".format(density)
        
        self.density = density
        self.mu      = mu
        self.bta     = bta
        
    def pdf(self, x):
        
        if self.density == Prior.GAUSSIAN:
            return np.sqrt(self.bta / (2. * np.pi)) * np.exp(-0.5 * np.power(x - self.mu, 2.))

        elif self.density == Prior.LOGNORMAL:
            return np.sqrt(self.bta / (2. * np.pi)) / x * np.exp(-0.5 * np.power(np.log(x / self.mu), 2.))

        
        
LMStepOutput = collections.namedtuple('LMStepOutput',
                                      ['yEst',
                                       'err',
                                       'WSS',
                                       'sigma2',
                                       'Objective'])

class LevenbergMarquardtReg: # frozen parameters
    
    def __init__(self, model_fn, lbda = .1, step_init = 1., min_displacement = 1E-5,
                 max_lbda = 1., step_mult_down = 0.8, step_mult_up = 1.2,
                 lbda_mult_up = 2., lbda_mult_down = 1.5,
                 check_every = 10, min_norm = 1E-5, max_iter = None):
        
        self.model_fn  = model_fn
        self.lbda      = lbda
        self.step_init = step_init
        self.min_displacement = min_displacement
        self.max_lbda  = max_lbda
        self.step_mult_down = step_mult_down
        self.step_mult_up   = step_mult_up 
        self.lbda_mult_up   = lbda_mult_up
        self.lbda_mult_down = lbda_mult_down
        self.check_every = check_every
        self.min_norm    = min_norm
        self.max_iter    = max_iter
        
        self.current_status = None
        
    def fit(self, X, y, theta_init, bounds = None, priors = None, weights = None):
        
        assert X.shape[0] == len(y), "Illegal input dimensions"
        
        self.nObs, self.nParams = X.shape
        
        self.X, self.y = X, y
        self.weights = np.ones(self.nObs) if weights is None else weights # not necessarily normalized
        self.lower, self.upper = self.__set_bounds__(bounds)
        self.priors = priors
                
        self.theta  = theta_init.copy()
        
        self.total_displacement = 0.
        self.step = self.step_init
        
        self.current_status = self.__get_optimization_status__(theta_init)
        print("Initial WSS: {}".format(self.current_status.WSS))
        
        nIter = 0
        while True:
            descent_direction = self.__find_descent_direction__()
            self.__move_to_new_theta__(descent_direction)
            nIter += 1
            if nIter % self.check_every == 0:
                norm_theta = np.linalg.norm(self.theta)
                if norm_theta == 0.:
                    raise Exception("Theta was set to 0")
                perc_displacement = self.total_displacement / norm_theta
                print("Check after {nIter} iterations: % displacement = {perc_displacement}, norm_theta = {norm_theta}" \
                      .format(nIter = nIter, perc_displacement = perc_displacement, norm_theta = norm_theta))
                
                if perc_displacement < self.min_norm:
                    break
                self.total_displacement = 0.

            if nIter == self.max_iter:
                break
            
    def predict(self, X):
        return self.model_fn(X, self.theta)

    def __find_descent_direction__(self):
        # descent direction solves a linear system Ax = b
        
        # Calculate descent direction from current theta
        JTWT = self.Jf_theta(self.X, self.theta)
        JTWT = np.dot(np.transpose(JTWT), np.sqrt(np.diag(self.weights)))
        A = np.dot(JTWT, np.transpose(JTWT)) # JT*WT*W*J
        b = np.dot(JTWT, self.current_status.err) # = - gradient of the objective function
        if self.priors is not None:
            A, b = self.__add_priors__(A, b)
            
        A += self.lbda * np.diag(np.diag(A)) # Marquardt
            
        success = False
        while True:
            if np.linalg.cond(A) < 1. / sys.float_info.epsilon:
                descent_direction = np.linalg.solve(A, b)
                success = True
                break
            if not self.__improve_conditioning__(A):
                break
            
        if not success:
            raise Exception("Could not calculate descent direction (singular matrix)")
            
        if np.dot(b, descent_direction) < -1E-10:
            raise Exception("Direction found is not a descent direction")
            
        return descent_direction
    
    def __add_priors__(self, A, b):
        
        A /= self.current_status.sigma2
        for j in range(len(self.priors)):
            pr = self.priors[j]
            if pr.density == Prior.GAUSSIAN:
                A[j][j] += pr.bta
                b[j]    -= pr.bta * (self.theta[j] - pr.mu)
            if pr.density == Prior.LOGNORMAL:
                log_theta_over_mu = np.log(self.theta[j] / pr.mu)
                A[j][j] += pr.bta / np.power(self.theta[j], 2.) #- (1. + (log_theta_over_mu - 1.) * pr.precision) / np.power(self.theta[i], 2.)
                b[j]    -= pr.bta * (log_theta_over_mu + 1. / pr.bta) / self.theta[j]

        return A, b

    def __move_to_new_theta__(self, descent_direction):
            
        norm_desc_dir = np.linalg.norm(descent_direction)
        descent_direction = descent_direction / norm_desc_dir
        
        self.status = self.__get_optimization_status__(self.theta)

        flg_theta_updated  = False
        while True:
            theta_new = self.theta + self.step * descent_direction
            if self.lower is not None:
                theta_new = np.clip(theta_new, self.lower, self.upper)
                
            new_status = self.__get_optimization_status__(theta_new)
            
            if new_status.Objective < self.current_status.Objective * (1. - 1E-5): # there has been a significant % decrease
                self.current_status = new_status
                self.theta = theta_new
                self.total_displacement += self.step * norm_desc_dir
                flg_theta_updated = True
                self.step *= self.step_mult_up
            else:
                if flg_theta_updated:
                    break
                self.step *= self.step_mult_down # try to decrease the step
            
            if self.step < self.min_displacement:
                break
        
        if not flg_theta_updated: # update lambda
            self.lbda = min(self.max_lbda, self.lbda * self.lbda_mult_up)
            
                
    def __get_optimization_status__(self, theta):
        
        yEst = self.model_fn(self.X, theta)
        err  = self.y - yEst
        WSS  = sum(self.weights * np.power(err, 2.))
        sigma2 = WSS # FIXME rivedere, va divisa per nObs per ottenere varianza stimata
        
        Objective = WSS
        if self.priors is not None:
            Objective = 0.5 * self.nObs * np.log(sigma2) + 0.5 * WSS / sigma2
            for j in range(len(self.priors)):
                pr = self.priors[j]
                if pr.density == Prior.GAUSSIAN:
                    Objective += 0.5 * pr.bta * np.power(theta[j] - pr.mu, 2.)
                elif pr.density == Prior.LOGNORMAL:
                    Objective += np.log(theta[j]) + 0.5 * pr.bta * np.power(np.log(theta[j] / pr.mu), 2.)
        
        return LMStepOutput(yEst      = yEst,
                            err       = err,
                            WSS       = WSS,
                            sigma2    = WSS,
                            Objective = Objective)
        
    def Jf_theta(self, X, theta, h = 1E-5):
        k = len(theta)
        
        Jf = []
        for i in range(len(theta)):
            Jf.append((self.model_fn(X, theta + h * np.eye(1, k, i)[0]) - self.model_fn(X, theta)) / h)
            
        return np.transpose(np.array(Jf))

    def __improve_conditioning__(self, A):
        
        flg_matrix_changed = False
        if max(abs(np.diag(A)) - 1.) > 1E-5:
            # Are there any zero rows in A? If so, put a 1. on their diagonal for those rows only.
            zero_rows = np.where(np.max(np.abs(A), axis = 1) < 1E-5)[0]
            if len(zero_rows) > 0:
                A.put([(A.shape[1] + 1) * i for i in zero_rows], 1.)
            else:
                # Last attempt: set all elements on the diagonal = 1.
                np.fill_diagonal(A, 1.)
                
            flg_matrix_changed = True

        return flg_matrix_changed
                
    def __set_bounds__(self, bounds):
        
        if bounds is None:
            lower = None
            upper = None
        else:
            lower, upper = [np.array(x) for x in zip(*bounds)]
            lower[lower == None] = -1E+30
            upper[upper == None] = +1E+30

        return lower, upper

if __name__ == '__main__':
    
    # Load training data
    X = train_data
    y = train_label
    X_train, y_train = X, y
    sorted_ixs = np.argsort(X[:, 0])
    theta_actual = np.array([15., .1, .4])
    
    # Define nonlinear model (Autoencoder) and declare LevenbergMarquardtReg class
    def f(X, theta):
        # return  X[:, 0])
        e = np.exp(1/(1+np.exp(np.dot(-theta[0], X) - theta[1])))
        d = np.dot(theta[0], e) + theta[2]
        f = np.linalg.norm(X - d, ord=2, keepdims=True)
        return f

    lr = LevenbergMarquardtReg(model_fn=f)
    lrReg = LevenbergMarquardtReg(model_fn=f)
    
    priors = [Prior(Prior.LOGNORMAL, mu = 20., bta = 10.), #Prior(Prior.GAUSSIAN, mu = 15., bta = 0.01),
              Prior(Prior.GAUSSIAN , mu = 1. , bta = 0.1),
              Prior(Prior.GAUSSIAN , mu = 1. , bta = 0.1)]
    
    # plot priors to get a qualitative idea
    fig, axs = plt.subplots(1, len(priors))
    for i in range(len(axs)):
        if priors[i].density == Prior.GAUSSIAN:
            sigma = 1. / np.sqrt(priors[i].bta)
            xx = np.linspace(priors[i].mu - 2. * sigma, priors[i].mu + 2 * sigma, 100)
        elif priors[i].density == Prior.LOGNORMAL:
            sigma = np.sqrt(1. / priors[i].bta)
            hlim = priors[i].mu * np.exp(sigma * (np.log(priors[i].mu) + 2 * sigma))
            xx = np.linspace(1E-5, hlim, 100.)
        axs[i].plot(xx, priors[i].pdf(xx))
        axs[i].set_title("parameter #{}".format(i), fontsize = 20)
    plt.suptitle("Prior distributions", fontsize = 25)
        
    theta_init = [np.mean(y_train), 0., 1.]
    print("theta_init = {}".format(theta_init))
    
    lr.fit(X_train, y_train, theta_init = theta_init)
    lrReg.fit(X_train, y_train, theta_init = theta_init, priors = priors)
    yEst_reg  = lr.predict(X)[sorted_ixs]
    
    # Display results
    expected_values = lr.__get_optimization_status__(theta_actual)
    print("\n*** RESULTS")
    print("Estimated theta:\n\t{} (reg)\n\t{} (non reg)".format(lrReg.theta, lr.theta))
    print("WSS for estimated theta:\n\t{} (reg)\n\t{} (non reg)".format(lrReg.current_status.WSS,
                                                                   lr   .current_status.WSS))
    print("WSS for real value of theta: {}".format(expected_values.WSS))
    
    plt.figure()
    plt.scatter(X[:, 0], y, color = 'black', marker = '.')
    plt.scatter(X_train[:, 0], y_train, color = 'red', marker = 'o', label= 'available observations')
    plt.legend()
    plt.title("Levenberg-Marquardt algorithm", fontsize=20)