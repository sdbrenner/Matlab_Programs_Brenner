function [binCents,avgB] = binAvg2(A,B,binSize,binRange,func)
% BIN_AVG Average data within bins that are defined in terms of a dependent
% variable.
%
%   avgA = binAvg(A,[],binsize) averages the data in 'A' within bins of
%   size binSize. Bin edges are defined by the start and end points of A,
%   and binSize:
%       binEdges = min(A):binSize:max(A);
%
%   [avgA,avgB] = binAvg(A,B,binSize) averages the data in both 'A' and 'B'
%   based on the dependent variable 'A'.  This allows for uneven spacing
%   between points, or uneven number of points present in each bin
%   (including the ability to have empty bins). 'A' must be a vector size
%   1xN or Nx1, and 'B' must be a vector or matrix with one dimension that
%   matches the length of A (i.e., B must be MxN or NxM), and cannot be a
%   scalar. If B is NxN, then I think averaging is performed along columns
%   of B (but results should be verified).
%
%   [avgA,avgB] = binAvg(...,binSize,binRange) allows for specification of
%   the start and end points of the averaging range for more control of the
%   bin edges:
%       binEdges = binRange(1):binSize:binRange(2);
%   
%   [avgA,avgB] = binAvg(...,func) allows for specification of the
%   averaging function, written in the form of an anonymous function.  By
%   default, func is '@nanmean', which is available from the statistics
%   toolbox.  If that toolbox isn't available, the same results can be
%   achieved by defining func as:
%       func = @(x) mean(x,'omitnan');
%
%   S.D.Brenner, 2019.


%% Error and input checking and organizing data

% 'A' must be a vector( 1xN or Nx1)
if ~isvector(A);  error('''A'' must be a vector'); end
% Find length of 'A'
N = length(A);

% rotate 'A' so it is 1xN 
rotA = false;       % keep track of original orientation
rotB = false;       % keep track of original orientation
if iscolumn(A)
    A = A';
    rotA = true;    % keep track of change in orientation  
end

% 'B' must be either a vector or a matrix with one dimension that matches 
% the length of 'A' (i.e., 'B' must be MxN or NxM), or it must be empty.
if ~isempty(B)
    if ~ismatrix(B); error('''B'' must be a vector or a matrix'); end
    [bM,bN] = size(B);
    if  bM ~= N && bN ~= N 
        error('The size of ''B'' must match the size of ''A'' in at least one dimension');
    end
    % Rotate 'B' so that it is MxN for consitency
    % ( averaging occurs along rows )
    if bM == N
        B = B';
        rotB = true;    % keep track of change in orientation    
    end
    % NOTE: If 'B' is a square matrix, this could end up causing problems
    % if averaging is meant to occur along the other dimension. Additional
    % support for such a use-case should be implemented here.
end


% if binRange is not specified, define it
if nargin < 4 || isempty(binRange); binRange = [min(A), max(A)]; end

% Check that binRange is a two-element vector
if length(binRange) ~= 2
    error('''binRange'' must be a 2-element vector'); 
end
binRange = sort(binRange);


% if func isn't specified, define it
if nargin < 5; func = @nanmean; end

%% Define Bins

% Bin size and edges
binEdges = binRange(1):binSize:binRange(end);

% For each data in 'A', find the corresponding bin number
bins = discretize(A, binEdges)';

% Identify NaN bins (which occur outside of binRange)
nanBins = isnan(bins);

% Reduce dimenion of 'A' and 'B' as appropriate
bins = bins(~nanBins);
A = A(~nanBins);
if ~isempty(B)
    B = B(:,~nanBins);
    [M,N] = size(B);  % overwrite previous value for N
else
    N = length(A);
end

%% Bin average (ignoring NaN's)

% For 'A' vector:
avgA = accumarray(bins,A,[],func,NaN).';

% For 'B': first check if B exists and is a vector or a matrix
if isempty(B)
    avgB = [];
elseif isvector(B)
    avgB = accumarray(bins,B,[],func,NaN).';
elseif ismatrix(B)
    [xx,yy] = ndgrid(1:M, bins );
    subs = [xx(:) yy(:)];
    rshpB = reshape( B, [1,M*N] );
    avgB = accumarray(subs,rshpB,[],func,NaN);
end

%% Re-categorize into defined bins

binCents = mean( [binEdges(1:end-1) ; binEdges(2:end) ]);
uBins = unique(bins);


evenB = NaN(M,length(binCents));
evenB(:,uBins) = avgB(:,uBins);
avgB = evenB;


%% Clean up data

% rotate 'avgA' and 'avgB' to match the original orientations of 'A' and 'B'
if rotA; avgA = avgA'; end
if rotB; avgB = avgB'; end

end