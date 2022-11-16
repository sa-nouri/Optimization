%% radiation treatment planning

close all; clear; clc;

%% 

run treatment_planning_data

cvx_begin
    variable b(n);
    minimize (sum(square_pos(Aother * b - Dother)))
    subject to
        b >= 0;
        b <= Bmax;
        Atarget * b >= Dtarget;
        %Atarget * b <= Dother;
cvx_end

%% plot

d_values = Atarget * b;
d_other = Aother * b;

figure(1)
histogram(d_values)
axis([0 2 0 60])
hold on
grid on
grid minor
plot([Dtarget Dtarget], [0 60], 'r')
hold off
title('Tumor voxel dose histogram')
xlabel('Dosage')

figure(2)
histogram(d_other)
axis([0 2 0 150])
hold on
plot([Dother Dother], [0 150], 'r')
grid on
grid minor
hold off
title('Other voxel dose histogram')
xlabel('Dosage')

%%