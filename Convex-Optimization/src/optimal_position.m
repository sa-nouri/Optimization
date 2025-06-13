%% Optimal political positioning

close all; clear; clc;

%% optimal political positioning

run opt_pol_pos_data

%% Counter example

plot_counterexample = true;
if plot_counterexample
    x = (-4:.01:4)';
    v_counterexample = -1;
    figure(1)
    plot(x, g(x + v_counterexample) + g( -x + v_counterexample))
    xlabel('x')
    ylabel('V')
    grid on
    grid minor
    print -deps opt_pol_pos_counterexample.eps
end

%% Approximation quality

plot_apx_quality = true;
if plot_apx_quality
    x = 0:.01:6;
    figure(2)
    plot(x, gapx(x), 'r', x, g(x), 'b')
    legend('gapprox', 'g', 'Location', 'Best')
    xlabel('z')
    grid on
    grid minor
    print -deps opt_pol_pos_logit_approx.eps;
    
    figure(3)
    plot(x, gapx(x) - g(x))
    ylabel('gapproxminusg')
    xlabel('z')
    grid on
    grid minor
    print -deps opt_pol_pos_logit_error.eps;
end

%% Compute optimal positions

cvx_begin quiet
    variable x(n)
    maximize(P' * gapx(W*x + v))
    subject to
        x >= l;
        x <= u;
        W*x + v >= 0;
cvx_end

% Vote before position optimization
Vi = P'*g(v)
Viapx = P'*gapx(v)

% Vote after position optimization
Vf = P' * g(W*x + v)
Vfapx = P' * gapx(W*x + v)

% Change in vote per constituency
delta = [P .* g(v), P .* g(W*x+v)]
