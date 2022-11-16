%% Amplifier designing

close all; clear; clc;

%% Computation

n = 4; 
Atot = 10000; 
alpha = [1e-5 1e-2 1e-2 1e-2]'; 
M = [0.1 5 10 10]'; 
Amax = [40 40 40 20]'; 

cvx_begin gp
    variables a(n) S(n) Sin
    prod(a) == Atot; 
    a <= Amax; 

    Nsquare(1) = a(1)^2 * alpha(1)^2;
    for i=2:n
        Nsquare(i) = a(i)^2*(Nsquare(i-1)+alpha(i)^2);
    end
   
    S(1) == a(1) * Sin;
    for i=2:n
        S(i) == a(i) * S(i-1);
    end
    
    S <= M; 
    maximize (S(n)/sqrt(Nsquare(n)))
cvx_end

opt_dyn_range = cvx_optval
opt_gains = a
signal_levels_and_limits = [S M]

%%