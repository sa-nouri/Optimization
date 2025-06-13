% robust_power_data.m
%
% Data file for the robust power problem

m = 20;  % Number of power stations
n = 30;  % Number of destination nodes

L = zeros(m, n);  % L[i, j] is length between station i and destination node j

% Lengths between nodes related to distance between i, j, but with some
% modification to keep solutions a bit more unique.
for ii = 1:m
  for jj = 1:n
    a = ii - .5;
    b = (jj - 1) / 2;
    mult = 2 + (ii - 1) / m;
    dist = 0;
    if a >= b
        dist = sqrt(a - b);
    else
        dist = sqrt(mult * (b - a));
    end
    L(ii, jj) = 1 + dist;
  end
end

% Transmission loss rate
alpha = .15;

% Power capacities on each power station
c = 5 * ones(m, 1);

% Usage level at node
u = 1.0 + linspace(0, 1, n);

