function xCents = binCenters(xEdges)
% BINCENTERS returns a vector of bin centers given a vector of bin edges
% (e.g. from 'histogram')
%
%   xCents = binCenters(xEdges)
%   xEdges must be an Nx1 vector, then xCents will be an (N-1)x1 vector
%
%   S.D.Brenner, 2019

%% Input parsing/error checking
% every good function should have this, but this is a lazy function so this
% will need to be added in after

xEdges = xEdges(:).';

%% Function execution
% really simple
xCents = mean( [xEdges(1:end-1) ; xEdges(2:end) ]);

end
