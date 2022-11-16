%% SDP solution for approximating non-convex least square problem

close all; clear; clc;
%%

randn('state', 0)
s_values = [0.5, 1, 2, 3];
m = 50;
n = 40;
A = randn(m,n);
xhat = sign(randn(n, 1));

for i = 1:length(s_values)
    s = s_values(i);
    b = A*xhat + s*randn(m, 1);
    
    cvx_begin sdp
        variable z(n, 1);
        variable Z(n, n) symmetric;
        minimize( trace(Z*A'*A) - 2*b'*A*z + b'*b )
        subject to
            [Z, z; z' 1] >= 0;
            diag(Z)==1;
    cvx_end
    fprintf(" s equals to : %f\n", s); 
    fprintf("The objective value for s equals is : %f\n", cvx_optval);
end

%% Compute f(a)

i = 1;      
s = s_values(i);
b = A*xhat + s*randn(m, 1);

cvx_begin
    variable x_ls(n, 1)
    minimize (pow_pos(norm(A*x_ls-b, 2), 2))
cvx_end

f_a = norm(A*sign(x_ls) -b);

%% compute f(b)

x_b = sign(z);
f_b = norm(A*sign(x_b) - b);

%% compute f(c)

[V, D] = eig(Z);
x_c = sign(V(:, length(Z)));
f_c = norm(A * x_c - b);

%% compute f(d)

mu = z;
Sigma = Z-(z*z');
rng('default')  % For reproducibility
R = mean(mvnrnd(mu, Sigma, 100),1)';
x_d = sign(R);
f_d = norm(A*x_d - b);

%%

P = [0, 0, 1/2, 1/4, 1/4, 0, 0;
    0, 0, 1/3, 0, 2/3, 0, 0;
    0, 0, 0, 0, 0, 1/3, 2/3;
    0, 0, 0, 0, 0, 1/2, 1/2;
    0, 0, 0, 0, 0, 3/4, 1/4;
    1/2, 1/2, 0, 0, 0, 0, 0;
    1/4, 3/4, 0, 0, 0, 0, 0];

yek = ones(length(P), 1);

cvx_begin
    variable p(length(P), 1)
    minimize( norm(P*p - p))
    subject to
        yek'*p==1;
cvx_end