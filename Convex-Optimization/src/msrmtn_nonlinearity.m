%% Measurment Non-linearity

run nonlin_meas_data

row = zeros(1, m);   row(1) = -1;    row(2) = 1;
col = zeros(1, m-1); col(1) = -1;

B = toeplitz(col, row);

cvx_begin quiet
    variables x(n) z(m);
    minimize(norm(z-A*x));
    subject to
        1/beta * B * y <= B * z;
        B * z <= 1/alpha * B * y;
cvx_end

fprintf('The estimated values of x are');
x
plot(z,y)
ylabel('y')
xlabel('z')
grid on
title('Maximum Likelihood estimation of \phi')

%%