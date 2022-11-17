%% Dodging Yogi Bear

close all; clear; clc;

%% Load the dataset

run hey_hey_booboo_data.m

%% Initial parameters

I = eye(2);
A = [I, delta*I; zeros(2,2), (1-alpha)*I];
B = [zeros(2); beta*I];
sigma_2 = 1/200;

X = zeros(size(V));
P_k = zeros(4);

for j = 0: (k-1)
    P_k = P_k + (A^j) * (B * B.') * (A^j)';
end

P_k = sigma_2 .* P_k;
P_k1 = P_k(1:2,1:2);


for i = 1:k
    X(i,:) = x_init_bear - ((1/norm(P_k1^0.5 * V(i, :).', 2)) ...
                        * threshold_05^0.5 * P_k1 * V(i, :).').';
end

%% Problem(c)

cvx_begin quiet
    variables u(k, 2) z(k+1, 4) x(k+1, 2)
    minimize(sum(pow_pos(u(:), 2)))
    
    subject to
        x(k, :) == x_safe;
        x(1, :) == x_init_human;
        z(1, :) == [ x_init_human, zeros(1, 2)];
        
        for i = 1: k
            z(i+1,:) == [A * z(i, :).' + B * u(i, :).'].';
            x(i+1, :) == ([eye(2), zeros(2)] * z(i,:).').';
            norm(u(i, :), 2) <= u_max;
            V(i, :) * ( X(i, :) - x(i+1, :) ).' >= 1;
        end
        
cvx_end

%% plot direction

figure(1)
plot(x(:, 1), x(:, 2), 'b')
title('Direction of human')
grid on
grid minor
hold on
plot(X(:, 1), X(:, 2), 'r')
legend("human's direction", "ellipsoid")

%% Problem (d)

v = zeros(size(V));

for i = 1: k
    v(i, :) = (x_init_bear - x(i, :))/norm(x_init_bear - x(i, :), 2);
end


cvx_begin quiet
    variables u(k, 2) z(k+1, 4) xc(k+1, 2)
    minimize(sum(pow_pos(u(:), 2)))
    
    subject to
        xc(k, :) == x_safe;
        xc(1, :) == x_init_human;
        z(1, :) == [ x_init_human, zeros(1, 2)];
        
        for i = 1: k
            z(i+1,:) == [A * z(i, :).' + B * u(i, :).'].';
            xc(i+1, :) == ([eye(2), zeros(2)] * z(i,:).').';
            norm(u(i, :), 2) <= u_max;
            v(i, :) * ( X(i, :) - xc(i+1, :) ).' >= 1;
        end
        
cvx_end


%% plot direction

figure(2)
plot(xc(:, 1), xc(:, 2), 'b')
title('Direction of human')
grid on
grid minor
hold on
plot(X(:, 1), X(:, 2), 'r')
plot(x(:, 1), x(:, 2), 'black')
legend("obtained V", "ellipsoid", "By given V")

%%