function [FX,FY] = m_gradient(F,lon,lat)
% M_GRADIENT Approximate gradient on lat,lon coordinates
%
%   [FX,FY] = m_gradient(F,lon,lat) calculates the gradient dF/dx, dF/dy
%   (where dx and dy are in [m]) for data definite on a lat-lon grid,
%   assuming a spherical earth of radius 6378.137km.
%
%   This function doesn't do any error checking and assumes "nice" lon, lat
%   inputs
%
%   Computed as:
%   (for lambda the longitude in radians and phi the latitude in radians)
%   d/dx = d(lambda)/dx * d/d(lambda);
%   d/dy = d(phi)/dy * d/d(phi);
%   where d/d(phi) = 1/R and d/d(lambda) = 1/(R*cos(phi))
%
%   S.D.Brenner, 2021

% Earth radius (same as is used in Rich's m_lldist in the M_Map package)
R = 6378.137 *1e3; % units of [m]

% Convert lat,lon to radians
lambda = deg2rad(lon);
phi = deg2rad(lat);

% Calculate transformations
DLambdaDx = 1./(R.*cos(phi));
DPhiDy = 1./R;

% Calculate and convert gradient
[DFDLambda,DFDPhi] = gradient(F,lambda,phi);
FX = DFDLambda.*DLambdaDx;
FY = DFDPhi.*DPhiDy;

end