%% 

close all; clear; clc;

%% Anlyzing of Flux Balance

% Load data
run fba_data

%% part(a): find maximum growth rate

cvx_begin quiet
    variable v(n)
    dual variable lambda
    maximize (v(n))
    subject to
        S*v == 0;
        v >= 0;
        lambda: v <= vmax;
cvx_end

Gstar = cvx_optval;
vopt = v;
[v vmax lambda]

%% part (b): find essential genes and synthetic lethals

G = zeros(n,n);

for i = 1:n
    for j = i:n
        cvx_begin quiet
            variable v(n)
            maximize (v(n))
            subject to
                S*v == 0;
                v >= 0;
                v <= vmax;
                v(i) == 0;
                v(j) == 0;
        cvx_end
        
        G(i,j) = cvx_optval;
        G(j,i) = cvx_optval;
    end
end

G < 0.2 * Gstar

%%
