%% DCP Rule check

close all; clear; clc;

%%
a = rand(1, 30);

cvx_begin
    variables x y u v z
    minimize( x + y + u + z) % Dummy Objective Function
    
    % Part(a)
    x == 0;
    y == 0;
    
    % Part(b)
    square_pos(square(x+y)) <= x -y 
    
    % Part(c)
    inv_pos(x) + inv_pos(y) <= 1;
    
    % Part(d)
    norm([v; u]) <= 3*x + y;
    max(x,1) <= v;
    max(y,2) <= u;
    
    % Part(e)
    x >= inv_pos(y);
    x >= 0;
    y >= 0;
    
    % Part(f)
    quad_over_lin(x + y, sqrt(y)) <= x - y + 5;
    
    % Part(g)
    quad_pos_over_lin(square(x), x) + quad_pos_over_lin(square(y), y) <= 1;
    pow_pos(x,3) + pow_pos(y,3) <= 1;
    
    % Part(h)
    x+z <= 1+geo_mean([x-quad_over_lin(z,y), y]);
    
cvx_end

%%