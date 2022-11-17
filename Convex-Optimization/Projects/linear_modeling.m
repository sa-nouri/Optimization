%% Robustness in linear modeling

close all; clear; clc;

%% Computation

run robust_linear_models_data

X_files = ["X_std", "X_outliers"];
Y_files = ["y_std", "y_outliers"];
Losses = ["sq", "abs", "norm"];

for i = 1: length(X_files)
    x = eval(X_files(i));
    for j = 1: length(Y_files)
        y = eval(Y_files(j));
        for k = 1: length(Losses)
            loss = Losses(k);
   
            cvx_begin quiet
                variable theta_opt(n)
                if strcmp(loss, "sq")
                   minimize(1/2 * (norm(x * theta_opt - y, 2)) )
                elseif strcmp(loss, "abs")
                    minimize(norm(x * theta_opt - y, 1))
                elseif strcmp(loss, "norm")
                    minimize(1./(max(norm(x,2), 1)) * norm(x * theta_opt - y, 2) )
                end
            cvx_end
                
            fprintf("The errors for %s, %s, and loss %s are \n", X_files(i), Y_files(j), Losses(k));
            errors = (theta_opt - theta_gen).'
            
            fprintf("The MABE for %s, %s, and loss %s is \n", X_files(i), Y_files(j), Losses(k));
            MABE = norm(errors, 2)/n
        end
    end
end
