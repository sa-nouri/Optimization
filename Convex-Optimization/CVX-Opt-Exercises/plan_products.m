%% Planning production with uncertain demand

close all; clear; clc;

run planning_data

%% Part(b) - Method(I)

cvx_begin quiet
    variables r(m) q(n)
    S = min(q* ones(1, K), D);
    maximize(p'* S * pi - c'*r)
    subject to
        r >= 0;
        q >= 0;
        r >= A*q;
cvx_end

rI = r;
disp('The optimal value of r')
disp(r);
disp('profit') 
disp(num2str(cvx_optval));
profitI = p.' * S - c.' * r;

%% Part(b) - Method(II)

cvx_begin quiet
    variables r(m) q(n, K)
    S = min(q, D);
    maximize(p'* S * pi - c'*r)
    subject to
        r >= 0;
        q >= 0;
        r*ones(1,K) >= A*q;
cvx_end

rII = r;
disp('The optimal value of r')
disp(rII);
disp('profit') 
disp(num2str(cvx_optval));
profitII = p.' * S - c.' * r;

