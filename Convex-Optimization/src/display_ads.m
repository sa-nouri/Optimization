%% Advertisement Displaying

close all; clear; clc;

%% 

run ad_disp_data

cvx_begin quiet
    variable N(n, T)
    s = pos(q-diag(Acontr' * N * Tcontr));
    maximize(R(:)' * N(:) - p' * s)
    subject to
        N >= 0;
        sum(N)' == I;
cvx_end

opt_s = s;
opt_net_profit = cvx_optval
opt_penalty = p' * opt_s
opt_revenue = opt_net_profit + opt_penalty

cvx_begin quiet
    variable N(n, T)
    maximize(R(:)' * N(:))
    subject to
        N >= 0;
        sum(N)' == I;
cvx_end

greedy_s = pos(q-diag(Acontr' * N * Tcontr));
greedy_net_profit = cvx_optval-p' * greedy_s
greedy_penalty = p' * greedy_s
greedy_revenue = cvx_optval

%%