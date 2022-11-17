## Numerical Optimization

Numerical optimization is defined as minimizing/maximizing an evaluation measure what is called “objective function” with several design constrains.

#### Gradient Based Methods

- [Gradinet Descent & Newton's method](GD_NM.py): Minimize the Rosenbrock function with steepest descent and Newton methods while using backtracking line search. Make comparison between the results.

- [Simulated Annealing](./simulated_annealing.py): It is the Simulated Annealing implementation, which is a stochastic global search optimization algorithm. This means that it makes use of randomness as part of the search process. This makes the algorithm appropriate for nonlinear objective functions where other local search algorithms do not operate well.

#### Conjugate Gradient Methods

- [Linear Conjugate Gradient Method](./conj_grad.py) and : It is the implementation of the conjugate gradient method which is an algorithm for the numerical solution of particular systems of linear equations, namely those whose matrix is positive-definite. The conjugate gradient method is often implemented as an iterative algorithm, applicable to sparse systems that are too large to be handled by a direct implementation or other direct methods such as the Cholesky decomposition. The [Cojugate Gradient Implementation](./conj_grad_linear_eig.py) is the conjugate gradient method for the following two different scenarios: one with eigen- values of A having uniform distribution and the other with eigenvalues of A being clustered around three different values.

- [Conjugate Gradient for AE](./conj_grad_ae.py): It is the implementation of the Fletcher-Reeves and Polak-Ribiere nonlinear conjugate gradient methods with back-tracking line search and also use the L-BFGS method to minimize the cost function of an auto- encoder network for the MNIST dataset considering.

#### Trust Region Methods

Trust-region method (TRM) first defines a region around the current best solution, in which a certain model (usually a quadratic model) can, to some extent, approximate the original objective function. TRM then take a step forward according to the model depicts within the region. Unlike the line search methods, TRM usually determines the step size before the improving direction (or at the same time).

- [Cauchy and Dog-Leg Methods](./cauchy_dogleg.py): It is the solution of the sub-problems using Cauchy point and Dog-Leg methods and comparison of their difference in performance.

- [Levenberg Marquardt](./levenberg_marquardt.py): It is the implementation of the Levenberg-Marquardt method to minimize the cost function of the auto-encoder network for the MNIST dataset.

#### Interior Point Methods

Interior point methods are a type of algorithm that are used in solving both linear and nonlinear convex optimization problems that contain inequalities as constraints.

- [Interior Piont Method SVM](./ipp_svm.py): It is goal is to optimize a soft margin linear SVM to classify the MNIST data using interior-point method.

#### Inexact Augmented Lagrangian Framework

- [iALM](./iALM/iALM.py): It is the implementation of inexact augmented Lagrangian method (iALM) for nonconvex problems with nonlinear constraints. It characterizes the total computational complexity of the method subject to a verifiable geometric condition, which is closely related to the Polyak-Lojasiewicz and Mangasarian-Fromovitz conditions. The paper of this implementation is available on [NIPS 2019 - iALM](https://proceedings.neurips.cc/paper/2019/file/866c7ee013c58f01fa153a8d32c9ed57-Paper.pdf)