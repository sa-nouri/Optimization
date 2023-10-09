# Convex Optimization

Convex optimization is a branch of mathematical optimization focused on minimizing convex functions over convex sets or maximizing concave functions over convex sets. Unlike many other optimization problems, certain types of convex optimization problems can be solved in polynomial time. 

This field has broad applications across various domains, including automatic control systems, estimation and signal processing, communications and networks, electronic circuit design, data analysis and modeling, finance, statistics (optimal experimental design), and structural optimization. Recent developments in computation and optimization algorithms have made convex programming more accessible and as understandable as linear programming.

The folder contains implementations of various problems, including:

### Basics

- **Objective Verification**: Implements simple objective functions to verify results against analytical solutions.

- **Check DCP**: Checks the understanding of various composition rules, convex analysis, and constraint reformulation rules.

- **Control Fuel Consumption**: Addresses the minimum fuel optimal control problem, where inputs are chosen to minimize total fuel consumption.

- **Boolean Linear Program**: Solves Boolean linear programs, where variables are constrained to 0 or 1.

- **Brute Force Constraint**: Provides a brute force solution to solving convex problems obtained by imposing subsets of constraints.

- **SDP for non-convex LS**: Focuses on non-convex least-squares approximation with binary constraints using semidefinite programming (SDP).

- **Data-BF-Constraints**: Contains data for the minimum number of constraints problem.

### Approximation and Fitting

- **Team Competition**: Models team competition outcomes as a convex optimization problem to estimate team abilities.

- **Team Data**: Provides data for team competition.

- **Measurement Non-linearity**: Estimates the maximum likelihood of an objective variable and its function, a problem involving infinite-dimensional maximum likelihood estimation.

- **Non-linear Measurement Data**: Contains necessary data for the measurement non-linearity problem.

- **Image Interpolation**: Focuses on image interpolation by minimizing roughness subject to interpolation conditions.

- **TV Image**: Defines matrices and data for TV image interpolation.

### Fitting

- **Numerical Filter Design**: Explores the design of a low-pass filter using convex optimization to meet filter specifications.

- **Fitting Sphere**: Minimizes the sum of square errors to fit a sphere to data.

- **Sphere Data**: Contains data for fitting a sphere.

- **Relational Function Fitting**: Fits a rational function to given data while constraining the denominator polynomial to be positive on an interval.

- **Optimal Operation of a Microgrid**: Optimizes the operation of a microgrid over time, including battery charge and locational marginal price (LMP) calculations.

- **Microgrid Data**: Provides data for optimal microgrid operation.

- **Amplifier Designing**: Optimizes the gains of amplifiers in a chain considering noise and overload effects.

- **Outlier Identification**: Identifies outliers using Löwner-John ellipsoids.

- **Anomaly Detection Data**: Contains data for outlier identification.

### Numerical Methods

- **Barrier Method**: Implements a barrier method for solving standard form linear programming (LP) and Newton's method for solving the centering problem.

- **Unconstrained Problem**: Solves unconstrained optimization problems using gradient descent and Newton's method.

- **Different Newton's Method**: Investigates variations of Newton's method to improve efficiency for large problems.

### Convex Optimization Problems

- **Disk Selection**: Chooses a set of disks to minimize an objective subject to constraints.

- **Disks Data**: Provides data for disk selection optimization.

- **Radiation Treatment**: Optimizes radiation delivery to minimize damage to healthy tissue while treating a tumor.

- **Treatment Data**: Contains data for radiation treatment optimization.

- **Flux Analyzing**: Focuses on flux balance analysis to compute bounds on reaction rates within a cell.

- **Flux Data**: Provides data for flux analyzing.

- **Advertisement Displaying**: Maximizes revenue by displaying ads optimally.

- **Ad Display Data**: Contains data for ad displaying optimization.

- **Optimal Position**: Finds optimal positions for political positioning.

- **Generator Dispatching**: Optimizes generator output power to meet electrical demand while considering dynamic constraints.

- **Generator Dispatch Data**: Provides data for generator dispatching.

### Convex Optimization Exercises

- **Optimizing Risk of Portfolio**: Addresses the portfolio risk-return trade-off problem.

- **Optimizing Risk of Uniform Problem**: Finds minimum-risk portfolios with the same expected return as the uniform portfolio.

- **Planning Productions**: Optimizes expected profit in different optimization problems with varying constraints.

- **Vehicle Lane Changing**: Plans a lane change for a vehicle.

- **Bicycle Optimal Trajectory**: Finds the optimal trajectory for a bicycle in various constraints.

### Projects

- **Linear Modeling Robustness**: Investigates different loss functions in linear predictor fitting.

- **Kernel SVM**: Implements the dual problem of kernel support vector machines.

- **Robust Powerline Provisioning**: Attaches powerlines for robustness against plant failures.

- **Dodging Yogi Bear**: Plans movement to avoid a bear while moving to a safe position.

- **Matrix Sensing**: Solves the matrix sensing problem using convex approximations.

This folder offers a wide range of practical problems and solutions in convex optimization. For detailed algorithmic explanations and implementation details, please refer to the respective code files within the directory. If you have any specific questions about any of these problems or algorithms, feel free to ask.
