%% Disk Selection

close all; clear; clc;

%% First Problem: Minimizing area

run disks_data

cvx_begin quiet
    variables C(n, 2) R(n)
    minimize(pi * sum(R.^2))
    subject to
        R >= 0;
        C(1:k, :) == Cgiven(1:k, :);
        R(1:k) == Rgiven(1:k);
        
        for i = 1 : size(Gindexes, 1)
            norm(C(Gindexes(i, 1), :) - C(Gindexes(i, 2), :)) ...
            <= (R(Gindexes(i, 1)) + R(Gindexes(i, 2)));
        end
        
cvx_end

fprintf('The optimal area: %3.4f\n', pi * cvx_optval);

%% Second Problem: Minimizing perimeter

cvx_begin quiet
    variables C(n, 2) R(n)
    minimize(2 * pi * sum(R))
    subject to
        R >= 0;
        C(1:k, :) == Cgiven(1:k, :);
        R(1:k) == Rgiven(1:k);
        for i = 1 : size(Gindexes, 1)
            norm(C(Gindexes(i,1), :) - C(Gindexes(i,2), :)) ...
            <= (R(Gindexes(i, 1)) + R(Gindexes(i, 2)));
        end
cvx_end
fprintf('optimal perimeter: %3.4f\n', 2 * pi * cvx_optval);

%% Plotting

% centers: n-by-2 matrix of circle centers
% radii: length-n vector of radii

centers = C;        radii = R;
figure(1)
viscircles(centers(1:k, :), radii(1:k), 'EdgeColor', 'r')
hold on
grid on
grid minor
viscircles(centers(k+1:n, :), radii(k+1:n), 'EdgeColor', 'b')
scatter(centers(1:k, 1), centers(1:k, 2), 100, 'ro', 'fill')
scatter(centers(k+1:n, 1), centers(k+1:n, 2), 100, 'bo', 'fill')

for i = 1 : size(Gindexes, 1)
    x1 = centers(Gindexes(i, 1),:); a = x1(1); b = x1(2);
    x2 = centers(Gindexes(i, 2),:); c = x2(1); d = x2(2);
    plot([a, c], [b, d], 'k-', 'LineWidth', 2);
end

xlim([-12.5, 12.5])
ylim([-12.5, 12.5])
pbaspect([1, 1, 1])
hold off
