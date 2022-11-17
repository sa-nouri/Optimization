## Convex Optimization

A branch of mathematical optimization known as "convex optimization" deals with the issue of minimizing convex functions over convex sets (or, equivalently, maximizing concave functions over convex sets). While mathematical optimization is typically NP-hard, certain types of convex optimization problems allow polynomial-time techniques.

Automatic control systems, estimation and signal processing, communications and networks, electronic circuit design, data analysis and modeling, finance, statistics (optimal experimental design), and structural optimization are just a few of the fields where convex optimization has applications. These fields have successfully applied the approximation concept. Convex programming is now practically as simple to understand as linear programming thanks to recent developments in computation and optimization algorithms.

The folder contains the following implementations:

#### Basics

- [Objective Verification](./Basics/obj_verification.m): It is about the implementations of simple objective functions to check the obtained resutls with the analytical solutions.

- [Check DCP](./Basics/dcp_check.m): This looks like a problem about ‘how to use CVX software’, or ‘tricks for using CVX’. But it really checks whether you understand the various composition rules, convex analysis, and constraint reformulation rules.

- [Control Fuel Consumption](./Basics/fuel_opt_control.m): The minimum fuel optimal control problem is to choose the inputs u(0), ..., u(N − 1) so as to minimize the total fuel consumed. This means that fuel use is proportional to the absolute value of the actuator signal, for actuator signals between −1 and 1; for larger actuator signals the marginal fuel eﬀiciency is half.

- [Boolean Linear Program](./Basics/boolean_linear_program.m): In a Boolean linear program, the variable x is constrained to have components equal to zero or one: minimize cT x subject to Ax ⪯ b xi ∈ {0,1}, i = 1,...,n
In general, such problems are very diﬀicult to solve, even though the feasible set is finite (containing at most 2n points). In a general method called relaxation, the constraint that xi be zero or one is replaced with the linear inequalities 0 ≤ xi ≤ 1.

- [Brute Force Constraint](./Basics/brute_force_subset_constraint.m): It is the brute force solution to solve all k convex problems obtained by choosing subsets of k constraints to impose, and selecting one with smallest objective value.

- [SDP for non-convex LS](./Basics/sdp_ls_approximator.m): It is about the non-convex least-squares approximation problem with binary constraints by deriving the dual of the problem in the form SDP.

- [Data-BF-Constraints](./Basics/satisfy_some_constraints_data.m): Data for satisfying minimum number of constraints problem


#### Approximation and Fittin

- [Team Competition](./Approximation/ml_team_competition.m): A set of n teams compete in a tournament. The problem is to find the maximum likelihood estimate of team abilities, given the outcomes, as a convex optimization problem. Then, we use the maximum likelihood estimate parameter found to predict the outcomes of next year’s tournament games, given in the matrix test, using yˆ(i) = sign(aˆj(i) − aˆk(i) ), and compare these predictions with the actual outcomes.

- [Team Data](./Approximation/team_data.m): This matrix gives the outcomes for a tournament in which each team plays each other team once.

- [Measurement Non-linearity](./Approximation/msrmtn_nonlinearity.m): The objective is to use convex optimization to find a maximum likelihood estimate of objective variable, as well as the function. This is an infinite-dimensional ML estimation problem, since one of the parameters we are estimating is a function.
Assumption: The measurements are numbered so that each output is sorted in nondecreasing order.

- [Non-linear Measurment Data](./Approximation/nonlin_meas_data.m): The necessary data for Measurement Non-linearity](msrmtn_nonlinearity.m) problem, which includes matrix A.

- [Image Interpolation](./Approximation/interpolate_image.m): The main focus is to interpolate the image by guessing the missing values. The reconstruction is found by minimizing a roughness measure subject to the interpolation con- ditions. One common roughness measure is the l2 variation (squared). Another method minimizes instead the total variation, which both models are convex optimization problem.

- [TV Image](./Approximation/tv_img_interp.m): This will define m, n, and matrices Uorig and Known. The matrix Known is m × n, with (i, j) entry
one if (i,j) ∈ K, and zero otherwise.

#### Fitting

- [Numerical Filter Design](./Fitting/design_filter.m): Consider the (symmetric, linear phase) finite impulse response (FIR) filter described by its fre- quency response. It is called a low-pass filter since low frequencies are allowed to pass, but frequencies above the cutoff frequency are attenuated. In this implementation we will explore the design of a low-pass filter using CVX to find the shortest length filter that satisfies the filter specification.

- [Fitting Sphere](./Fitting/fit_sphere.m): The implementation is to minimize the sum of square error to fit a sphere to data.

- [Sphere Data](./Fitting/sphere_fit_data.m): Contains the necessary data for the implementation of fitting sphere problem.

- [Relational Function fitting](./Fitting/fit_relational_func.m): In this implementation we fit a rational function p(t)/q(t) to given data, while constraining the denominator polynomial to be positive on the interval [α, β]. The optimization variables are the numerator and denominator coefficients ai, bi. The interpolation points ti ∈ [α,β], and desired function values yi,i = 1,...,k, are given.

- [Optimal Operation of a Microgrid](./Fitting/op_microgrid.m): We consider a small electrical microgrid that consists of a
photovoltaic (PV) array, a storage device (battery), a load, and a connection to an external grid. We will optimize the operation of the microgrid over one day, in 15 minute increments, so all powers, and the battery charge. Then, we want to find the locational marginal price (LMP). The LMPs can be used as a system for payments among the load, the PV array, the battery, and the grid.

- [Microgrid Data](./Fitting/microgrid_data.m): The necessary data for the implementation fo the optimal operation of a microgrid problem.

- [Amplifier Designing](./Fitting/design_amplifier.m): The system consists of n amplifiers connected in a chain. The variables that we will optimize over are the gains of the amplifiers. The two concerned effects are noise generated by the amplifiers, and amplifier overload. 

- [Outlier Identification](./Fitting/Identify_outlier.m): This implementation is about an outlier identification technique using L ̈owner-John ellipsoids. Roughly speaking, we think of an outlier as a point that is far away from most of the points.

- [Anomaly Detection Data](./Fitting/ellip_anomaly_data.m): It is the anomaly detection data for identifying outlier by the technique implemented in [Outlier Identification](./Fitting/Identify_outlier.m)

#### Numerical Methods

- [Barrier Method](./Numerical/barrier_method.m): It is the implementation a barrier method for solving the standard form LP as well as Newton’s method for solving the centering problem. It also includes LP solver with strictly feasible starting point for solving the standard form of LP.

- [Unconstrained Problem](./Numerical/unconstrained_problem.m): It is about implementation of the gradient descend method and Newton's method for solving a unconstrained problem with using reasonable choices for the backtracking parameters, and a stopping criterion.

- [Different Newton's Method](./Numerical/newtons_method.m): The cost of Newton’s method is dominated by the cost of evaluating the Hessian ∇2f(x) and the cost of solving the Newton system. For large problems, it is sometimes useful to replace the Hessian by a positive definite approximation that makes it easier to form and solve for the search step.
    1. Re-using the Hessian. We evaluate and factor the Hessian only every N iterations, where N > 1, and use the search step ∆x = −H−1∇f(x), where H is the last Hessian evaluated.
    2. Diagonal approximation. We replace the Hessian by its diagonal, so we only have to evaluate the n second derivatives, and computing the search step is very easy.

#### Convex Optimization Problems

- [Disk Selection](./CVX-Opt-Problems/select_disk.m): The goal is to choose a set of n disks (i.e., specify their centers and radii) to minimize an objective subject to some constraints.

- [Disks Data](./CVX-Opt-Problems/disks_data.m): The necessary data for solving the optimization problem of disk selection.

- [Radiation Treatment](./CVX-Opt-Problems/perform_radiation.m): In radiation treatment, radiation is delivered to a patient, with the goal of killing or damaging the cells in a tumor, while carrying out minimal damage to other tissue. The radiation is delivered in beams, each of which has a known pattern; the level of each beam can be adjusted. In most cases multiple beams are delivered at the same time, in one ’shot’, with the treatment organized as a sequence of ’shots’. A minimum radiation dose should be administered to each tumor voxel.

- [Treatment Data](./CVX-Opt-Problems/treatment_planning_data.m): Here we have split the matrix A into Atarget, which contains the rows corresponding to the target voxels, and Aother, which contains the rows corresponding to other voxels. It is mandatory for solving [Radiation Treatment](./CVX-Opt-Problems/perform_radiation.m) problem.

- [Flux Analyzing](./CVX-Opt-Problems/analyze_flux.m): Flux balance analysis is based on a very simple model of the reactions going on in a cell, keeping track only of the gross rate of consumption and production of various chemical species within the cell. Based on the known stoichiometry of the reactions, and known upper bounds on some of the reaction rates, we can compute bounds on the other reaction rates, or cell growth, for example.

- [Flux Data](./CVX-Opt-Problems/fba_data.m): The data for the implementation of the [Flux Analyzing](./CVX-Opt-Problems/analyze_flux.m).

- [Advertisement Displaying](./CVX-Opt-Problems/disp_ads.m): The problem is about to maximize revenue for dispalying ads, we would simply display the ad with the highest revenue per impression, and no other, in each display period. [Ad Display Data](./CVX-Opt-Problems/ad_disp_data.m) is the data for implementing this optimization problem.

- [Optimal Position](./CVX-Opt-Problems/optimal_position.m): A political constituency is a group of voters with similar views on a set of political issues. The electorate (i.e., the set of voters in some election) is partitioned (by a political analyst) into K constituencies, with (nonnegative) populations. A candidate in the election has an initial or prior position on each of n issues, but is willing to consider (presumably small) deviations from her prior positions in order to maximize the total number of votes she will receive. This implementation is to fiind the optimal positions for the partisan political positioning problem with data given in [Optimal Position Data](./CVX-Opt-Problems/opt_pol_pos_data.m).

- [Generator Dispatching](./CVX-Opt-Problems/generator_dispatch.m): In the generator dispatch problem, we schedule the electrical output power of a set of generators over some time interval, to minimize the total cost of generation while exactly meeting the electrical demand. One challenge in this problem is that the generators have dynamic constraints, which couple their output powers over time. For example, every generator has a maximum rate at which its power can be increased or decreased. It is the solution of the generator dispatch problem with the data given in [Generator Dispatch Data](./CVX-Opt-Problems/gen_dispatch_data.m), which gives (fake, but not unreasonable) demand data for 2 days, at 15 minute intervals. This file includes code to plot the demand, optimal generator powers, and prices.

#### Convex Optimization Exercises

- [Optimizing Risk of Portfolio](./CVX-Opt-Exercises/optmize_risk_pf.m): It is about the portfolio risk-return trade-off problem. The main focus is to obtain the optimal portfolios and compute the expected return and standard deviation. 
This portfolio maximizes the expected return, subject to the probability of a loss being no more than 5%. The monte carlo simulation comes to solve this problem in this case.

- [Optimizing Risk of Uniform Problem](./CVX-Opt-Exercises/optmize_risk_pf_uniform.m): The implementation is about to find minimum-risk portfolios with the same expected return as the uniform portfolio, with risk measured by portfolio return variance, and different portfolio constraints, with data that can be found in the file [Portfolio Data](./CVX-Opt-Exercises/simple_portfolio_data.m).

- [Planning Productions](./CVX-Opt-Exercises/): We want to explore two different optimization problems that arise in choosing the variables. The objective is to maximize the expected profit. Then, the goal is to carry out the methods from part (a) on the problem instance with numerical data given in [Planning Data](./CVX-Opt-Exercises/planning_data.m). This file will define A, D, K, c, m, n, p and pi. The K columns of D are the possible demand vectors.

- [Vehicle lane changing](./CVX-Opt-Exercises/change_lane.m): A vehicle is traveling down a highway with two lanes. The goal of this implementation is to plan a lane change.

- [Bicycle Optimal Trajectory](./CVX-Opt-Exercises/opt_trajectory.m): It is about finding the optimal trajectory of the bicycle in different constraints.

#### Projects

- [Linear Modeling Robustness](./Projects/linear_modeling.m): In this simulation, we investigate the choice of losses in a problem of fitting a linear predictor to given data. We consider three losses, each giving different robustness properties: the squared error, the absolute error, and the normalized error. 
The data for this problem is available in [Robustness Linear Model Data](./Projects/robust_linear_models_data.m). There are two data matrices X and two target vectors y in the file, Xstd, Xoutliers, ystd, and youtliers. The pair Xstd, ystd corresponds to data generated via the well-specified model (1), with w a mean-zero vector. The matrix Xoutliers has its first 10 rows corrupted by large noise, and similarly, the vector youtliers has its first 10 entries corrupted.

- [Kenrnel SVM](./Projects/kernel_svm.m): Using the kernel K(x,z) = (1+xTz)^6 and objective f(t) = max{0, 1−t}, the file is the implementation of the dual problem of kernel support vector machine for the problem data in [Kernel SVM DATA](./Projects/kernel_svm_data.m).

- [Robust powerline provisioning](./Projects/powerline_provising.m): We wish to attach powerlines to minimize the cost of provisioning the lines while protecting against plant failures. We represent these lines as a bipartite graph between the m plants and n destination nodes. Then, it is unrealistic to assume that all powerstations and powerlines will remain up, so you would like to build robustness into your system. We would like to guarantee m − k reliability: the system will provide the desired level of power uj at each destination node even if k of the power plants fail. Using the [Robust Power Data](./Projects/robust_power_data.m), the file is the implementation of this problem in CVX.

- [Dodging Yogi Bear](./Projects/find_actions.m): A human and a bear (Yogi) begin at specific positions. The bear wanders aimlessly, but if we and the bear get too close, it will steal our picnic basket and eat it; we wish to avoid this outcome (at least with high probability). Our goal is to move to a safe position while avoiding the bear. Using the file [Boo Boo Data](./Projects/hey_hey_booboo_data.m), we want to implement our solution to the optimizatin problem.

- [Matrix Sensing](./Projects/matrix_sensing.m): In the matrix sensing problem, the observations come from in a particular form. This is an instance of a composite optimization problem, which we can (approximately) minimize by sequentially minimizing convex approximations. Using the data in the file [Matrix Data](./Projects/matrix_sco_data.m), which defines a U matrix, V matrix, and b vector, the goal is to implement the dual form of the problem.. Have the  procedure iterate until the change between iterations satisfies a constraint.
