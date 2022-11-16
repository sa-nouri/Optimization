function ellipse_draw(A,b)

% make a unit ball    
teta = 0:pi/50:2*pi;
unit_ball = [cos(teta);sin(teta)];

%inverse of unit ball under an affine function
ellips = A\(unit_ball-b);

plot(ellips(1,:),ellips(2,:))

end
