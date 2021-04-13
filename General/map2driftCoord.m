function [xd,ld] = map2driftCoord(lon,lat,theta,lon0,lat0,lonRange,latRange)
% MAP2DRIFTCOORD Converts lon, lat to local "drift" coordinates
%
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
%
%   [xd,ld] = map2driftCoord(lon,lat,theta,lon0,lat0,lonRange,latRange)
%   allows specification of map range being considered.  Points in lon, lat
%   outside of this range will return NaN values for drift coordinates. The
%   inputs lonRange,latRange must both be monotonic vectors of length 2.
%   
%
%   Notes: 
%   (1) This function makes use of the M_Map toolbox made by Rich
%   Pawlowicz, and uses that to do the "heavy lifting" in terms of
%   calculation.  This function was devloped based on M_Map v1.4. (M_Map is
%   availble free at: https://www.eoas.ubc.ca/~rich/map.html ) 
%   (2) Calculation is done by first converting to a set of UTM coordinates
%   in the map area specified by lonRange,latRange (or by the min,max range
%   of the input lon,lat vectors).  The UTM conversion is most accurate at
%   the center of the map area, and errors increase further away (this
%   means that the same set of lat,lon points can translate into slight
%   different drift coordintes if the mapping area differs).
%
%   S.D.Brenner, 2018-2019

%% Parse inputs

% Modify longitude convension:
% lons should measure from -180 to 180 degrees (instead of 0 to 360) 
lon = mod(lon+180,360)-180;


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

% If lonRange, latRange are not specified determine appropriate limits
if nargin < 6 || isempty(lonRange) || isempty(latRange)
    latRange = [min([lat(:);lat0]),max([lat(:);lat0])];
    lonRange = [min([lon(:);lon0]),max([lon(:);lon0])];
elseif numel(lonRange) ~= 2 || numel(latRange) ~= 2
    error('lonRange and latRange must be monotonic vectors of length 2');
end

% The orgin point must be within the range
if ~( latRange(1) <= lat0 && latRange(2) >= lat0 ) ||...
   ~( lonRange(1) <= lon0 && lonRange(2) >= lon0 )
    error('Origin point must be within the specified mapping area');
end



%% Convert to UTM coordinates
tempFig = figure; m_proj('UTM','lat',latRange,'lon',lonRange); close(tempFig);
% tf = figure; m_proj('UTM','lat',[71,75],'lon',[-150,-140]); close(tf);
[X,Y] = m_ll2xy(lon,lat);

% Define origin point in UTM coordinates
[X0,Y0] = m_ll2xy(lon0,lat0);

%% Rotate and convert to across and along-drift directions

% Define direction of along-drift vector:        
th = mod(450-theta,360); %[deg] % converstion to cartesian coordinate system

% Convert from UTM to LD and XD components
ld =  (X-X0)*cosd(th) + (Y-Y0)*sind(th) ;
xd = -(X-X0)*sind(th) + (Y-Y0)*cosd(th) ;

% Convert from [m] to [km]
ld = ld/1000;
xd = xd/1000;