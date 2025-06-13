%% Numerical Filter Design

close all; clear; clc;

%% Computation coefficients

K = 500;
wp = pi/3; wc = .4 * pi; alpha = 0.0316;
w = linspace(0, pi, K);
wi = max(find(w <= wp)); wo = min(find(w >= wc));

for N = 1:50
    k = [0:1:N]';
    C = cos(k*w)';
    cvx_begin
        variable a(N+1)

        C(1:wi, :) * a <= 1.12;
        C(1:wi, :) * a >= 0.89;
        cos(wp * linspace(0, N, N+1)) * a >= 0.89;
        
        C(wo:K, :) * a <= alpha;
        C(wo:K, :) * a >= -alpha;
        cos(wc * linspace(0, N, N+1)) * a <= alpha;
    cvx_end
    if (strcmp(cvx_status,'Solved') == 1)
        break;
    end
end

%% Filter Computation and plotting

H = a'*cos(k*w);
plot(w, H)
set(gca, 'YTick', [-alpha 0 alpha .89 1 1.12])
axis([0 pi -alpha 1.12])
grid on
hold on
plot([0 wp wp], [.89 .89 -alpha], ':')
plot([wc wc pi], [1.12 alpha alpha], ':')
xlabel('\omega')
ylabel('H(\omega)')
title("Designed Low pass filter")
hold off
%%
