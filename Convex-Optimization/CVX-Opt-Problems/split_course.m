%% Theory-applications split in a course

close all; clear; clc;

%% theory-applications split in a course

a = 2; b = 3;
theta = [ 0.05; 0.1; 0.3];  alpha = [-0.1; 0.8; -0.3];
beta = [ 1.4; -0.3; 0.7];

n = 20;     m = 3; 

for plan = 1:m+1
    cvx_begin quiet
        variable T(n)
        expressions s(m, n+1) obj

        for i = 1:n
            s(:, i+1) = (1-theta).*s(:, i) + theta .* (alpha * T(i) + beta * (1 - T(i)));
        end

        if plan == 4
            obj = min(s(:, n+1));
        else
            obj = s(plan, n+1);
        end
        
        maximize obj
        subject to
            T >= 0;
            T <= 1;
            T(1:b) == 1;
            cumsum(1-T(b+1:n)) <= a * cumsum(T(b+1:n));
    cvx_end

    subplot(4, 1, plan);
    
    plot(1:n, T, 'k', ...
    1:n, s(1, 2:n+1), 'r', ...
    1:n, s(2, 2:n+1), 'g', ...
    1:n, s(3, 2:n+1), 'b');

    title(sprintf('Plan %d', plan));
    grid on
    grid minor
    axis([1 n -0.5 1.5]);
    xlabel('Lectures')
    fprintf('Plan %d: %f %f %f\n', plan, s(1, n+1), s(2, n+1), s(3, n+1));
end