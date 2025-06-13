%% Robust powerline provisioning

close all; clear; clc;

%% Part (a)

run robust_power_data.m

cvx_begin quiet
    variables A(m, n) P_tr(m, n) P_ar(m, n)
    minimize(trace(L.' * A))
    subject to
        for i=1:m
            p_tr = 0;
            for j=1:n
                P_tr(i, j) >= P_ar(i, j) + (alpha) * L(i, j) * quad_over_lin(P_ar(i, j), A(i, j));
                p_tr = p_tr + P_tr(i,j);
            end
            p_tr <= c(i);
        end
        for i=1:n
            p_ar = 0;
            for j=1:m
                p_ar = p_ar + P_ar(j,i);
            end
            p_ar >= u(i);
        end
        P_tr >= 0;
        P_ar >= 0;
        A(:) >= 0;
cvx_end

fprintf('\nThe total cost of solution is: %f \n', cvx_optval)

fprintf('\nThe total elements of A, which are greater than 0.001, is: %i \n',...
    length(find(abs(A) > 0.001) ))

%% Part (b)
% Implementing robust optimization

K = [1, 2, 3, 4];

for k=1:length(K)
    
    cvx_begin quiet
        variables A(m, n) P_tr(m, n) P_ar(m, n)
        minimize(trace(L.' * A))
        subject to
            for i=1:m
                p_tr = 0;
                for j=1:n
                    P_tr(i, j) >= P_ar(i, j) + (alpha) * L(i, j) * quad_over_lin(P_ar(i, j), A(i, j));
                    p_tr = p_tr + P_tr(i,j);
                end
                p_tr <= c(i);
            end
            for i=1:n
                sum_smallest(P_ar(:, i), m-k) >= u(i);
            end
            P_tr >= 0;
            P_ar >= 0;
    cvx_end
    
    fprintf('\nThe results for k equals: %i \n', k)
    fprintf('\nThe total cost of solution is: %f \n', cvx_optval)
    fprintf('\nThe total elements of A, which are greater than 0.001, is: %i \n',...
        length(find(abs(A) > 0.001) ))

end

%%
