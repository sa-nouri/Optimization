%% Team competition prediction

close all; clear; clc;

% Estimation

run team_data

A = sparse(1:m, train(:,1), train(:,3), m, n) + ...
  sparse(1:m, train(:,2), -train(:,3), m, n);

cvx_begin quiet
    variable a_est(n)
    minimize(-sum(log_normcdf(A*a_est/sigma)))
    subject to
        a_est >= 0;
        a_est <= 1;
cvx_end

fprintf('The estimated values of a are\n');
disp(a_est');

% Prediction

A_test = sparse(1:m_test, test(:,1), 1, m_test, n)+ ...
     sparse(1:m_test, test(:,2), -1, m_test, n);

residuals = sign(A_test * a_est);

p_maximum_likelihood = 1 - length(find(residuals - test(:, 3)))/m_test;
train_same_out = 1 - length(find(train(:, 3) - test(:, 3)))/m_test;

fprintf('The maximum likelihood prediction accuracy is %f \n', p_maximum_likelihood);
fprintf('The same outcome of train is %f \n', train_same_out);