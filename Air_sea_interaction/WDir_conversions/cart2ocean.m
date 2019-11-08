function [DIR,R] = cart2ocean(U,V,convention,units)
%CART2OCEAN transforms cartesian to a version of polar coordinates with the
% convention for direction consistent with commonly used oceanographic choices.
%
%   [DIR,R] = cart2ocean(U,V) transforms corresponding elements of data
%   stored in Cartesian coordinates U,V to polar coordinates (angle DIR
%   and magnitude R).  The arrays U and V must be the same size (or
%   either can be scalar). DIR is returned in degrees using the
%   'oceanographic' convention for direction, which is the direction the
%   vector is pointing TO on a compass (0 = North, 90 = East, 180 = South,
%   270 = West).  The units of R match the units of U and V.
%
%   [DIR,R] = cart2ocean(U,V,convention) allows for the choice of direction
%   conventions. 'convention' can be one of 'nautical', 'oceanographic', or
%   'polar'. 
%
%   The 'oceanographic' convention is as described above (compass
%   direction TO), and is commonly used for descriting ocean currents.
%
%   The 'nautical' convention is the direction that the vector is pointing
%   FROM on a compass (0 = North, 90 = East, 180 = South, 270 = West).
%   This is the convention commonly used in describing wind or wave
%   directions.
%
%   The 'polar' convention is that used in mathematics, so is the direction
%   that the vector is pointing TO, measured clockwise from the x-axis 
%   (0 = East, 90 = North, 180 = West, 270 = South).
%
%   [DIR,R] = cart2ocean(...,units) allows for the choice of output units
%   between radians and degrees. 'units' must be one of 'rad' or 'deg'.
%
% S. Brenner, 2018

%% Set defaults
if nargin < 4 || isempty(units); units = 'deg'; end
if nargin < 3 || isempty(convention); convention = 'oceanographic'; end


%% Convert from u,v to polar cooridinates
[theta,R] = cart2pol(U,V); 
theta = rad2deg(theta);

%% Adjust direction values based on choice of convention
switch convention
    case 'polar'
        DIR = mod(theta,360);
    case 'nautical'
        DIR = mod(270-theta,360);
    case 'oceanographic'
        DIR = mod(90-theta,360);
    otherwise
        error('The input variable ''convention'' must be one of ''nautical'', ''oceanographic'', or ''polar''');
end

%% Adjust output units
switch units
    case 'deg'
    case 'rad'
        DIR = deg2rad(DIR);
    otherwise
        warning('The input variable ''units'' must be one of ''rad'' or ''deg''.  Using ''deg''.');
end

end

