%%
close all; clear; clc;

%% x1 + x2
cvx_begin
    variable x(1, 2)
    minimize( sum(x)) % part(a)
    subject to
        (2*x(1) + x(2)) >= 1;
        (x(1) + 3*x(2)) >= 1;
         x >= 0;
cvx_end
print(x)

%% -(x1 + x2)
cvx_begin
    variable x(1, 2)
    minimize( -sum(x)) % part(b)
    subject to
        (2*x(1) + x(2)) >= 1;
        (x(1) + 3*x(2)) >= 1;
         x >= 0;
cvx_end
print(x)

%% x1
cvx_begin
    variable x(1, 2)
    minimize(x(1)) % part(c)
    subject to
        (2*x(1) + x(2)) >= 1;
        (x(1) + 3*x(2)) >= 1;
         x >= 0;
cvx_end
print(x)
 
%% max(x)

cvx_begin
variable x(1, 2)
    minimize(max(x)) % part(d)
    subject to
        (2*x(1) + x(2)) >= 1;
        (x(1) + 3*x(2)) >= 1;
         x >= 0;
cvx_end
print(x)

%% x(1)^2 + 9*(x(2)^2)

cvx_begin
variable x(1, 2)
    minimize(x(1)^2 + 9*(x(2)^2)) % part(e)
    subject to
        (2*x(1) + x(2)) >= 1;
        (x(1) + 3*x(2)) >= 1;
         x >= 0;
cvx_end
print(x)

%%