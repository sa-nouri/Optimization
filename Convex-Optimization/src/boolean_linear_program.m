%% ‌Boolean Linear Program

close all; clear; clc;

% Generate problem data

rand('state', 0);
n = 100;  m = 300;
A = rand(m, n);      b=A*ones(n, 1)/2;    c = -rand(n, 1);

% solve LP relaxation
cvx_begin
    variable x(n)
    minimize (c'*x)
    subject to
        A*x <= b;
        x> = 0;
        x< = 1;
cvx_end

x_rlx = x;
L = cvx_optval;

% Examining the thresholds
thresholds = 0:0.01:1;  max_violations = zeros(length(thresholds), 1);
objectives = zeros(length(thresholds), 1);

for i=1:length(thresholds)
    x_hat = (x_rlx >= thresholds(i));
    max_violations(i) = max(A*x_hat - b);
    objectives(i) = c'*x_hat;
end

% least upper bound and its threshold
i_feasible = find(max_violations <= 0);
U = min(objectives(i_feasible));

t_m = min(i_feasible);
min_threshold = thresholds(t_m);

% plot objective and max violation versus threshold
subplot(2, 1, 1)
plot(thresholds(1: t_m-1), max_violations(1: t_m-1), 'r', thresholds(t_m:end), max_violations(t_m:end), 'b');
xlabel('threshold');
ylabel('max violation');
grid on

subplot(2,1,2)
hold on; plot(thresholds, L*ones(size(thresholds)), 'k');
plot(thresholds(1:t_m-1), objectives(1:t_m-1), 'r', thresholds(t_m:end), objectives(t_m:end), 'b');
xlabel('threshold');
ylabel('objective');
grid on
%%