function [lon,lat] = EASEgrid2map(r,s)
% [lon,lat] = EASEgrid2map(r,s)
%
% EASE-Grid Projection Coordinate tranformation
% Additional details: 
% http://nsidc.org/ease/clone-ease-grid-projection-gt
%
% Consistent with the definitions used for Polar Pathfinder sea-ice
% velocity (https://nsidc.org/data/nsidc-0116)

% Define geometry constants (from NSIDC)
R = 6371.228;   % Radius of the Earth
C = 25.067525;  % Nominal cell size
r0 = 180;       % Map origin column 
s0 = 180;       % Map origin row 

r = r - r0;
s = s - s0;

lambda = atan2(r,s);
phi = pi/2 - 2*asin( C*sqrt(r.^2+s.^2)/(2*R) );

lon = rad2deg(lambda);
lat = rad2deg(phi);

end