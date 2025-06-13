%% Interpolate Image by guessing the missing values

close all; clear; clc;
%%

run tv_img_interp

cvx_begin quiet
    variable Ul2(m, n);
    Ul2(Known) == Uorig(Known);
    Ux = Ul2(1:end, 2:end) - Ul2(1:end, 1:end-1);
    Uy = Ul2(2:end, 1:end) - Ul2(1:end-1, 1:end);
    minimize(norm([Ux(:); Uy(:)], 2));
cvx_end

cvx_begin quiet
    variable Utv(m, n);
    Utv(Known) == Uorig(Known);
    Ux = Utv(1:end,  2:end) - Utv(1:end, 1:end-1);
    Uy = Utv(2:end, 1:end) - Utv(1:end-1, 1:end);
    minimize(norm([Ux(:); Uy(:)], 1)); 
cvx_end 

%%
% Graph everything.
figure(1); cla;
colormap gray;

subplot(221);
imagesc(Uorig)
title('Original image');
axis image;

subplot(222);
imagesc(Known.*Uorig + 256-150*Known);
title('Obscured image');
axis image;

subplot(223);
imagesc(Ul2);
title('l_2 reconstructed image');
axis image;

subplot(224);
imagesc(Utv);
title('Total variation reconstructed image');
axis image;