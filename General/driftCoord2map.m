function [lon,lat] = driftCoord2map(xd,ld,theta);
% driftCoord2map Converts local "drift" coordinates to lat, lon
%	[lon,lat] = driftCoord2map(xd,ld) converts input arrays of "along-
%   drift" (ld) and "across-drift" (xd) coordinates to latitude, longitude
%	arrays using a coordinate rotation angle of 280 degrees
%
%   [lon,lat] = driftCoord2map(xd,ld,theta) allows specification of the
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

% Convert from [km] to [m]
ld = ld*1000;
xd = xd*1000;

% Define origin point in UTM coordinates
tf = figure; m_proj('UTM','lat',[71,75],'lon',[-150,-140]); close(tf);
[X0,Y0] = m_ll2xy(lon0,lat0);

% Define direction of along-drift vector        
th = mod(450-theta,360); %[deg] % conversion to cartesian coordinate system

% Convert from LD and XD to UTM components
X = X0 + ( ld*cosd(th) - xd*sind(th) );
Y = Y0 + ( ld*sind(th) + xd*cosd(th) );


% Convert from UTM coordinates to lat, lon
m_proj('UTM','lat',[71,75],'lon',[-150,-140]);
[lon,lat] = m_xy2ll(X,Y);


