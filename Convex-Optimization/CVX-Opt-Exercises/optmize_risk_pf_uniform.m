%% Risk-return trade-off in portfolio optimization

close all; clear; clc;

run simple_portfolio_data

%% part(a) -1: No(additional) constraint

cvx_begin
    variable x(n)
    minimize(quad_form(x, S))
%     minimize(x.' * S * x)
    subject to
        sum(x) == 1;
        pbar' * x == x_unif' * pbar;
cvx_end

x_without_constraint = x;

%% part(a) -2: Long-only

cvx_begin
    variable x(n)
    minimize(quad_form(x, S))
%     minmize(x.' * S * x)
    subject to
        sum(x) == 1;
        x >= 0;
        pbar' * x == x_unif' * pbar;
cvx_end

x_long_only = x;

%% Limit on total short position

cvx_begin
    variable x(n)
    minimize(quad_form(x, S))
%     minmize(x.' * S * x)
    subject to
        sum(x) == 1;
        sum(pos(-x)) <= 0.5;
        pbar' * x == x_unif' * pbar;
cvx_end

x_limit_position = x;

%% Optimial Risks in the obtained x-values

uniform = sqrt(quad_form(x_unif, S));
no_constraint_opt = sqrt(quad_form(x_without_constraint, S));
long_only_opt = sqrt(quad_form(x_long_only, S));
lim_short_position =  sqrt(quad_form(x_limit_position, S));

disp('The No (additional) constraint sd:');
disp(no_constraint_opt)

disp('The long only sd:')
disp(long_only_opt)

disp('The Limit on short position sd:')
disp(lim_short_position)

disp('The x_unif sd:')
disp(uniform)

%% Optimal risk-return trade-off curves

novals = 100;   mu_values = logspace(-1, 4, novals);
r_long_only = [];    r_lim_short_position = [];
sd_long_only = [];   sd_lim_short_position = [];

i = 1;
while i < novals
    mu = mu_values(i);
    
    cvx_begin
        variable x_long(n)
        maximize(pbar' *x_long - mu * quad_form(x_long, S))
        subject to
            sum(x_long)==1;
            x_long>=0;
    cvx_end
    
    r_long_only = [r_long_only, pbar' * x_long];
    sd_long_only = [sd_long_only, sqrt(x_long' * S * x_long) ];
    
    cvx_begin
        variables x_lim_pos(n)
        maximize(pbar' * x_lim_pos - mu * quad_form(x_lim_pos, S))
        subject to
            sum(x_lim_pos)==1;
            sum(pos(-x_lim_pos)) <= 0.5;
    cvx_end
    
    r_lim_short_position = [r_lim_short_position, pbar'*x_lim_pos];
    sd_lim_short_position = [sd_lim_short_position, sqrt(x_lim_pos'*S*x_lim_pos)];
    i = i +1;
end

%% plot figures

figure(1)
plot(sd_long_only, r_long_only, 'b');
hold on;
plot(sd_lim_short_position, r_lim_short_position, 'r');
grid on
grid minor
xlabel('standard deviation of portfolio return')
ylabel('mean return')
legend('Long only', 'Limit Short Position')
