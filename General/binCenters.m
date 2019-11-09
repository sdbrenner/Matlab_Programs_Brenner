function binCents = binCenters(binEdges)
% BINCENTERS finds bin centers
%
%   binCents = binCenters(binEdges) returns a vector of bin centers given a
%   vector of bin edges (e.g. from 'histogram'). binEdges must be a vector.
%   If binEdges is size 1xN or Nx1, then binCents is size 1x(N-1) or
%   (N-1)x1 respectively.
%
%   S.D.Brenner, 2019

%% Input parsing/error checking

% Check type
if ~isvector(binEdges)
    error('Input must be a vector');
end

% Ensure xEdges is a row vector
xEdgesRow = binEdges(:).';

    

%% Function execution
% ( this is really simple )

% Find bin centers
binCents = mean( [xEdgesRow(1:end-1) ; xEdgesRow(2:end) ]);

% Rotate to the correct orientation
if iscolumn(binEdges)
    binCents = binCents(:);
end


end
