%% Optimal Trajectory of the Bicycle

close all; clear; clc;

%% Solve problem

x_initial = [-5; 0]; x_final = [6; 1]; v_initial = [2; 0];
T = 12;     h = 0.1;    n = T/h;

cvx_begin
    variables x(2, n) v(2, n) a(2, n)
    minimize(h * sum(pow_pos(a(:), 2)))
    subject to
        for i = 1 : n -1
            x(:, i+1) == x(:, i) + h/2 * (v(:, i) + v(:, i+1));
            v(:, i+1) == v(:, i) + h*a(:, i);
        end
        x(:, 1) == x_initial;
        x(:, n) == x_final;
        v(:, 1) == v_initial;
cvx_end

x_prev = x;

%% Plot optimal trajectory

figure(1)
plot(x(1, :), x(2, :))
grid on
grid minor
title('Optimal trajectory of bicycle')
xlabel('x')
ylabel('y')

disp("The optimal cost of trajectory bicycle is")
disp(cvx_optval)

%% Solve Roundabout problem and plot problem

% x_prev = x;

c = zeros(size(x_prev));

for i = 1:n
    c(:, i) = x_prev(:, i)/norm(x_prev(:, i));
end

cvx_begin
    variables x(2, n) v(2, n) a(2, n)
    minimize(h * sum(pow_pos(a(:), 2)))
    subject to
        for i = 1 : n -1
            x(:, i+1) == x(:, i) + h/2 * (v(:, i) + v(:, i+1));
            v(:, i+1) == v(:, i) + h*a(:, i);
        end
        
        x(:, 1) == x_initial;
        x(:, n) == x_final;
        v(:, 1) == v_initial
        
        for i = 1:n
            c(:, i).' * x(:, i) >= 1
        end 
cvx_end

figure(2)
plot(x(1, :), x(2, :))
grid on
grid minor
title('Optimal trajectory of bicycle')
xlabel('x')
ylabel('y')

sprintf("The optimal cost of trajectory bicycle is \n")
disp(cvx_optval)

%% Trajectory converges

% x_prev = x;

cvx_opt_value = 0;
epsilon = 1e-2;

while(1)
    
    c = zeros(size(x_prev));
    for i = 1:n
        c(:, i) = x_prev(:, i)/norm(x_prev(:, i));
    end
    
    cvx_begin
        variables x(2, n) v(2, n) a(2, n)
        minimize(h * sum(pow_pos(a(:), 2)))
        subject to
            for i = 1 : n -1
                x(:, i+1) == x(:, i) + h/2 * (v(:, i) + v(:, i+1));
                v(:, i+1) == v(:, i) + h*a(:, i);
            end

            x(:, 1) == x_initial;
            x(:, n) == x_final;
            v(:, 1) == v_initial

            for i = 1:n
                c(:, i).' * x(:, i) >= 1
            end 
    cvx_end
    x_prev = x;
    
    if abs(cvx_optval - cvx_opt_value) < epsilon
        disp("Optimal value converged")
        break
    end
    cvx_opt_value = cvx_optval;
end

figure(2)
plot(x(1, :), x(2, :))
grid on
grid minor
title('Optimal ctrajectory covergence of bicycle')
xlabel('x')
ylabel('y')

disp("The optimal cost of trajectory convergence bicycle is ")
disp(cvx_optval)

%%