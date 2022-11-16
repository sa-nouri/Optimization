%% Optimal operation of a microgrid

close all; clear; clc;

run microgrid_data 

%% Solving Optimization

cvx_begin quiet
    variables p_batt(N) p_buy(N) p_sell(N) p_grid(N) q(N)
    variables payments income
    dual variable v
    minimize(payments - income)
    subject to
    % Grid constraints
    payments == (.25) * p_buy' * R_buy
    income == (.25) * p_sell' * R_sell
    p_grid == p_buy - p_sell;
    v : p_ld == p_grid + p_batt + p_pv;
    p_buy >= 0;
    p_sell >= 0;
    % Battery constraints
    q >= 0;
    q <= Q;
    p_batt <= D;
    p_batt >= -C;
    q(1) == q(N) - (.25)*p_batt(N);
    for i = 2: N
        q(i) == q(i-1) - (.25)*p_batt(i-1)
    end
cvx_end

fprintf('Optimal Value: $%.2f\n', cvx_optval);

%% Plot things

pp_start = 34;
p_start = 48;
p_end = 72;
pp_end = 86;

% Plot vertical lines where prices change and horizontal lines
% where the limits of each value are, for interpretability

figure()
plot(p_grid)
title('Grid Power (kW)');
title('Power (kW)');
xlabel('Interval');
line([pp_start, pp_start], [-35, 35], 'Color', 'black', 'LineStyle', '--');
line([p_start, p_start], [-35, 35], 'Color', 'black', 'LineStyle', '--');
line([p_end, p_end], [-35, 35], 'Color', 'black', 'LineStyle', '--');
line([pp_end, pp_end], [-35, 35], 'Color', 'black', 'LineStyle', '--');
line([1,N], [0, 0], 'color', 'black');
ylim([-35, 35])
figure()
plot(p_batt)
title('Battery Power (kW)');
ylabel('Power (kW)');
xlabel('Interval');
line([pp_start, pp_start], [-C, D], 'Color', 'black', 'LineStyle', '--');
line([p_start, p_start], [-C, D], 'Color', 'black', 'LineStyle', '--');
line([p_end, p_end], [-C, D], 'Color', 'black', 'LineStyle', '--');
line([pp_end, pp_end], [-C, D], 'Color', 'black', 'LineStyle', '--');
line([1,N], [D, D], 'color', 'black', 'linestyle', '--');
line([1,N], [-C, -C], 'color', 'black', 'linestyle', '--');
line([1,N], [0, 0], 'color', 'black');
ylim([-C-2, D+2])
figure()

plot(q)
title('Battery Charge (kWh)');
title('Energy (kWh)');
xlabel('Interval');
line([pp_start, pp_start], [0, Q], 'Color', 'black', 'LineStyle', '--');
line([p_start, p_start], [0, Q], 'Color', 'black', 'LineStyle', '--');
line([p_end, p_end], [0, Q], 'Color', 'black', 'LineStyle', '--');
line([pp_end, pp_end], [0, Q], 'Color', 'black', 'LineStyle', '--');
line([1,N], [Q, Q], 'color', 'black', 'linestyle', '--');
line([1,N], [0, 0], 'color', 'black');
ylim([0, Q+2])

%% Getting LMP

if v(1) < 0
    v = -v;
end

LMP = 4*v;

figure()
hold on
plot(R_buy, 'linestyle', '--')
plot(R_sell, 'linestyle', '--')
plot(LMP)
legend('Buy Price', 'Sell Price', 'LMP');
title('LMP vs Nominal Prices');
xlabel('Interval');
ylabel('Price ($/kWh)');
grid on

%% Calculating Payments

load_cost = p_ld' * v
batt_cost = -p_batt' * v
PV_cost = -p_p' * v
grid_costs = p_grid' * v
net_cost = -load_cost - batt_cost - PV_cost + grid_costs

%%
