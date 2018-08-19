function [U,V] = ocean2cart(DIR,R,convention,units)
%OCEAN2CART transforms oceanographic polar to cartesian coordinates where U
%is a is in the northward direction and V is in the eastward direction
%   [U,V] = ocean2cart(DIR,R) transforms corresponding elements of data
%   stored in polar coordinates (angle DIR, radius R) to Cartesian
%   coordinates X,Y.  The arrays DIR and R must the same size (or
%   either can be scalar).  DIR is in degrees using the 'oceanographic' 
%   convention for direction, which is the direction the vector is pointing 
%   TO on a compass (0 = North, 90 = East, 180 = South, 270 = West).
%
%   [U,V] = ocean2cart(DIR,R,convention) can be used if DIR is specified
%   using a convention other than oceanographic.  'convention' can be one of 
%   'nautical', 'oceanographic', or 'polar'.
%   The 'oceanographic' convention is as described above.  The 'nautical'
%   convention is the direction that the vector is pointing FROM on a
%   compass (0 = North, 90 = East, 180 = South, 270 = West).  The 'polar'
%   convention is that used in mathematics, so is the direction that the
%   vector is pointing TO (0 = East, 90 = North, 180 = West, 270 = South).
%
%   [U,V] = ocean2cart(...,units) can be used if DIR is specified in a set
%   of units other than degrees. 'units' must be one of 'rad' or 'deg'.
%
% S. Brenner, 2018



transforms corresponding elements of data
%   stored in Cartesian coordinates X,Y to polar coordinates (angle DIR
%   and magnitude R).  The arrays u and v must be the same size (or
%   either can be scalar). DIR is returned in degrees using the
%   'oceanographic' convention for direction, which is the direction the
%   vector is pointing TO on a compass (0 = North, 90 = East, 180 = South,
%   270 = West).
%



%% Set defaults 
% assuming input in in oceanograpic convention and units are degrees
if nargin < 4; units = 'deg'; end
if nargin < 3; convention = 'oceanographic'; end

%% Convert directions to a polar cooridinate system
switch convention
    case 'polar'
        theta = DIR;
    case 'nautical'
        theta = mod(270-DIR,360);
    case 'oceanographic'
        theta = mod(90-DIR,360);
    otherwise
        error('''convention'' must be one of ''nautical'', ''oceanographic'', or ''polar''');
end

%% Convert angles to radians
switch units
    case 'deg'
        theta = deg2rad(theta);
    case 'rad'
    otherwise
        warning('''units'' must be one of ''rad'' or ''deg''.  Assuming ''deg''.');
        theta = deg2rad(theta);
end

%% Convert to cartesian coordinates
[U,V] = pol2cart(theta,R);


end