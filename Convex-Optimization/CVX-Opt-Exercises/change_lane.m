%% Smoothest possible lane change

close all; clear; clc;

%% Define Parameters and Solve the problem

T = 30;     Tstart = 15;    Tend = 20;  time = 1:T;
Smin = 25; Smax = 35;
L = 3.7;

cvx_begin
    variables x(T+1) y(T+1)
    acceleration = [x(3:T+1); y(3:T+1)] - 2*[x(2:T);y(2:T)] + [x(1:T-1); y(1:T-1)];
%     velocity = norms([x(2:T+1); y(2:T+1)] - [x(1:T); y(1:T)], 2);
%     selected_velocity = velocity(Tstart:Tend);
    minimize(sum(pow_pos(acceleration(:), 2)))
    subject to
        y(1:Tstart + 1) == 0;
        y(Tend + 1:end) == L; 
        x(1) == 0;
        y(1) == 0;
        y <= L;
        y >= 0;
        x >= 0;
        x(2:end) - x(1:end-1) >= zeros(T, 1);
%         velocity <= ones(T, 1) * Smax;
        [x(2:Tstart), y(2:Tstart)] - [x(1:Tstart-1), y(1:Tstart-1)] >=  ones(Tstart-1, 2) * Smin;
        [x(Tend:end), y(Tend:end)] - [x(Tend-1:end-1), y(Tend-1:end-1)] >= ones(length(x)-Tend +1, 2) * Smin;
        [x(2:end), y(2:end)] - [x(1:end-1), y(1:end-1)] <= ones(length(x)-1, 2) * Smax;
cvx_end

%% plot velocity

figure(1)
plot(time, velocity)
grid on
grid minor
xlabel('Time [s]')
ylabel('Velocity [m/s]')

%% plot position

figure(2)
plot(x, y, 'b-*')
grid on
grid minor
xlabel('x')
title('Position')
ylabel('y')

%%

x = (1:31).';
y = (1:31).';
p = [x, y];

p(2:Tstart + 1, 1) - p(1:Tstart, 1) >= ones(Tstart, 1) * Smin;
p(Tend:end, 1) - p(Tend-1:end-1, 1) >= ones(12, 1) * Smin;