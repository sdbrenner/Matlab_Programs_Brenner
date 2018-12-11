function h = trackRibbon(x,y,z,C)
% TRACKRIBBON 
%
%   trackRibbon(x,y,z,C)
%   x and y are both 1xN, z is Mx1, and C is MxN.
%
%   h = trackRibbon(...) returns the function handle

%% Error checking
% Obviously this still needs to be added.  Mostly this should check input
% sizes and types.


%% Ensure data are arranged correctly
x = x(:)';
y = y(:)';
z = z(:);

%% Create variable matrices

[M,N] = size(C);

X = repmat(x,[M,1]);
Y = repmat(y,[M,1]);
Z = repmat(z,[1,N]);

%% Plot

h = surf(X,Y,-Z,C);
shading interp;

end