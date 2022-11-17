%% Matrix sensing and sequential convex optimization

close all; clear; clc;

run matrix_sco_data

%% Implement Procedure

num_initilizations = 10;
epsilon = 1e-4;
x_initial = randn(n_1, num_initilizations);
y_initial = randn(n_2, num_initilizations);

c_vals = {};

for idx = 1: num_initilizations
   
    x_vals = [];  y_vals =[];
    x = x_initial(:, idx);    y = y_initial(:, idx);
    
    k = 1;
    f_vals = [];

    while(1)
        
        for i = 1:m
            u = U(i, :);
            v = V(i, :);
            linear_mapping(i, 1) = trace(v.' * u * x * y.');
        end
        
        cvx_begin quiet
            variables delta_x(n_1) delta_y(n_2)
            minimize(norm(b - linear_mapping - ...
                    diag(V * y) * U * delta_x - diag(U * x) * V * delta_y, 1))
        cvx_end
        if (norm(delta_x, 2)^2 + norm(delta_y, 2)^2) <= epsilon^2
            break;
        end

        f_vals = [f_vals; cvx_optval];
        x_vals = [x_vals; x];
        y_vals = [y_vals; y];

        x = x + delta_x;
        y = y + delta_y;
    end
    
    disp(idx)
    c_vals{idx} = f_vals;
    
end



%% Plot Error values

index = 10;
figure(1)
semilogy(c_vals{index})
hold on
grid on
grid minor
xlabel('iterations')
ylabel('error')
title('f(x,y)')
hold off

%%