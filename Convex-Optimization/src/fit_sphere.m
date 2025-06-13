%% Fitting a sphere to data

close all; clear; clc;

%% Computation coefficients

run sphere_fit_data.m

A = [2*U.', ones(1,length(U)).'];
b = vecnorm(U, 2).';

cvx_begin
    variable x(3,1)
    cvx_precision best
    minimize(norm(A*x - b))
cvx_end

%% Center and radius of circle

r = sqrt(x(3) - norm(x(1:2)));
x_c = x(1:2); 

ang = 0:0.01:2*pi; 
xp = r * cos(ang);
yp = r * sin(ang);

%% Plot

figure()
scatter(U(1,:), U(2,:))
hold on
plot(x_c(1) + xp, x_c(2) + yp);
grid on;
xlabel('x')
ylabel('y')
title('Fitting sphere to data')
%%
