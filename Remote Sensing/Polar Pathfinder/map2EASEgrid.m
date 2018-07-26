function [r,s] = map2EASEgrid(lon,lat)
% [r,s] = map2EASEgrid(lon,lat)
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

lambda = deg2rad(lon);
phi =    deg2rad(lat);

r = 2*R/C * sin(lambda) .* sin(pi/4 - phi/2) + r0;
s = 2*R/C * cos(lambda) .* sin(pi/4 - phi/2) + s0;

end