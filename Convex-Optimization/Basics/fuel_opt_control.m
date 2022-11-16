%% Fuel optimal controlw

close all; clear; clc;

%%

N = 30;     state_dimension = 3;

A = [ -1 0.4 0.8; 1 0 0 ; 0 1 0];
b = [ 1; 0; 0.3];
x_des = [ 7; 2; -6];

x0 = zeros(state_dimension, 1); % initial state

cvx_begin
    variable x(state_dimension, N+1)
    variable u(1, N)
    minimize (sum(max(abs(u), 2*abs(u)-1)))
    subject to
        x(:, 2:N+1) == A*x(:, 1:N)+ b*u; % main constrain, Linear Dynamical Systems
        x(:, 1) == x0;      % initial state
        x(:, N+1) == x_des; % target stat
cvx_end

% Plot 
stairs(0: N-1, u, 'linewidth', 1.5)
title('Actuater Signal u(t)')
xlabel('time')
ylabel('u')
grid on