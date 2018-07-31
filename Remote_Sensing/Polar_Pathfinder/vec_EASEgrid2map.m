function [LON,LAT,UL,VL] = vec_EASEgrid2map(X,Y,U,V)
%
% [LON,LAT,UL,VL] = vec_EASEgrid2map(X,Y,U,V)
% 
[LON,LAT] = EASEgrid2map(X,Y);

[Th,S] = cart2pol(U,V);
Th = rad2deg(Th);
Th_off = theta_offset(X,Y);
Th_mod = mod( Th_off - Th, 360);

[UL,VL] = pol2cart( deg2rad(Th_mod) , S );
end


function th_off = theta_offset(X,Y)
%
% th_off = theta_offset(X,Y)
%
% Returns the direction offset to convert from EASE-grid theta to mapping
% theta.

[X0,Y0] = map2EASEgrid(0,90);

dX = X - X0;
dY = Y - Y0;
alpha = atan2d( dY,dX);
th_off = alpha-90;
th_off = mod(th_off,360);
end