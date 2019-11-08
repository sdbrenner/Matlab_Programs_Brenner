function [h,coords] = trackRibbon(x,y,z,C)
% TRACKRIBBON plots along a three-dimension "ribbon"
%
%   trackRibbon(x,y,z,C) plots values from the matrix C as colour along a
%   three-dimensional ribbon defined by x, y, and z, where x and y are both
%   1xN, z is Mx1, and C is MxN.
%
%   h = trackRibbon(...) returns the function handle
%   [h,coords] = trackRibbon(...) also returns a data structure, coords,
%   containing the gridded X,Y,Z coordinates along the track ribbon.

%% Error checking
% Obviously this still needs to be added.  Mostly this should check input
% sizes and types.


%% Ensure data are arranged correctly
x = x(:)';
y = y(:)';
z = z(:);

%% Create gridded variable matrices

[M,N] = size(C);

X = repmat(x,[M,1]);
Y = repmat(y,[M,1]);
Z = repmat(z,[1,N]);

%% Plot

hold on;
h = surf(X,Y,Z,C);
shading interp;

%% Save coordinates in data structure:

coords.X = X;
coords.Y = Y;
coords.Z = Z;

end