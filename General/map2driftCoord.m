function [xd,ld] = map2driftCoord(lon,lat,theta,lon0,lat0);
% MAP2DRIFTCOORD Converts lon, lat to local "drift" coordinates
%	[xd,ld] = map2driftCoord(lon,lat) converts longitude, latitute array
%	inputs and to "along-drift" (ld) and "across-drift" (xd) coordinate
%	arrays (measured in km) using a coordinate rotation angle of 280
%	degrees
%
%   [xd,ld] = map2driftCoord(lon,lat,theta) allows specification of the
%   rotation angle, theta
%
%   [xd,ld] = map2driftCoord(lon,lat,theta,lon0,lat0) allows specification
%   of the origin point


% If no theta is specified, use a default value:
if nargin < 3 || any(isempty(theta)) || any(isnan(theta))
    theta =  280;
end
% If no origin is specified, use a default value
if nargin < 5 || length(lon0) ~= 1 || length(lat0) ~= 1
    lon0 = -145.484;
    lat0 = 73.1426;
    
    lon0 = -145.1815264372435;
    lat0 = 73.124582298727731;
end


% Convert to UTM coordinates
mapLatRange = [min([lat(:);lat0]),max([lat(:);lat0])];
mapLonRange = [min([lon(:);lon0]),max([lon(:);lon0])];
tf = figure; m_proj('UTM','lat',mapLatRange,'lon',mapLonRange); close(tf);
% tf = figure; m_proj('UTM','lat',[71,75],'lon',[-150,-140]); close(tf);
[X,Y] = m_ll2xy(lon,lat);

% Define origin point in UTM coordinates
[X0,Y0] = m_ll2xy(lon0,lat0);

% Define direction of along-drift vector:        
th = mod(450-theta,360); %[deg] % converstion to cartesian coordinate system

% Convert from UTM to LD and XD components
ld =  (X-X0)*cosd(th) + (Y-Y0)*sind(th) ;
xd = -(X-X0)*sind(th) + (Y-Y0)*cosd(th) ;

% Convert from [m] to [km]
ld = ld/1000;
xd = xd/1000;