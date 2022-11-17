%% Kernel Support Vector Machines

close all; clear; clc;

%% Load and define parameters

run kernel_svm_data

epsilon = 1e-6;
K = (1 + X * X.').^6;
G = K + eye(size(K)) * epsilon;

rmin = -1.5;    rmax = 1.5;     n = length(G);
Lambda = [1e-3, 1e-2, 1, 10];
r = linspace(rmin, rmax, n)';
x_in = [r, r];
[r1p, r2p] = meshgrid(x_in, x_in);
data_gen = [r1p(:) r2p(:)];

%% Computation

lambda = Lambda(4);

cvx_begin
    cvx_precision best
    variable alpha_opt(length(G), 1)
%     maximize( -ones(length(G), 1).' * max(1-alpha_opt, 0) - (1/(2*lambda)) * quad_form(y.*alpha_opt, G) )
    maximize( -(ones(length(G), 1).' * alpha_opt) - (1/(2*lambda)) * quad_form(y.*alpha_opt, G))
    subject to
        alpha_opt <= 0;
        alpha_opt >= -1;
cvx_end

v = -1/lambda * diag(y) * alpha_opt;

kernel = (1+X*data_gen.').^6;
phi_theta = reshape(kernel'*v, size(r1p));

%% Data distribution

plot_scatter(X, y)
grid on
title('Scatter plot of data')

%% Plot Contour

figure(2)
contour(r1p, r2p, phi_theta);
title('Contour plot of {\phi_\theta} for {\lambda} = 10')
figure(3)
contourf(r1p, r2p, phi_theta, 'ShowText','on');
title('Contour plot of {\phi_\theta} for {\lambda} = 10')

%% Scatter plot function

function scatter_figure = plot_scatter(X, Y)
    oneindices = find(Y == 1);
    minusoneindices = find(Y == -1);
    Xones = X(oneindices, :);
    Xminusones = X(minusoneindices, :);
    scatter_figure = figure;
    scatter(Xones(:,1), Xones(:,2));
    hold on;
    scatter(Xminusones(:,1), Xminusones(:,2));
    hold on;
end