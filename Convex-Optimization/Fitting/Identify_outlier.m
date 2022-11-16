%% Outlier Identification Technique

close all; clear; clc;

%% ellipsoidal peeling

ellip_anomaly_data;
volumes = [];
removed = [];
figure(1); hold all;

for i = 1:6

    cvx_begin
        variable A(2,2) symmetric
        variable b(2)
        dual variable v
        maximize (det_rootn(A))
        subject to
            v : norms(A*X + b*ones(1,size(X,2))) <= 1;
    cvx_end
    ellipse_draw(A,b);

    [vm idx] = max(v);
    removed = [ removed X(:,idx) ];
    X(:,idx) = [];
    volumes = [ volumes 1/det(A) ];
end

plot(X(1,:), X(2,:), 'bx'); % normal points
plot(removed(1,:), removed(2,:), 'ro'); % outliers

figure(2);
semilogy(volumes)
xlabel('number of removed points');
ylabel('ellipsoid volume');
set(gca,'XTick', 1:length(volumes));

