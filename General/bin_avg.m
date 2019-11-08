function [avgA,avgB] = bin_avg(A,B,binsize,data_range,func)
% BIN_AVG Average data within bins.
%   [avg_a,avg_b] = bin_avg(a,b,binsize) averages the data in bins of size
%   'binsize' in the dependent variable 'a'.
%
%   [avg_a,avg_b] = bin_avg(a,b,binsize,data_range) allows for
%   specification of the start and end points of the averaging range
%   
%   [avg_a,avg_b] = bin_avg(a,b,binsize,data_range,func) allows for
%   specification of the averaging function, written in the form of an
%   anonymous function.  By default, func is '@nanmean'

%% Error and input checking and organizing data

% 'A' must be a vector( 1xN or Nx1)
if ~isvector(A);  error('''A'' must be a vector'); end
% rotate 'a' to row vector:
N = length(A);

% 'B' must be either a vector or a matrix with one dimension that matches 
% the length of 'A' (i.e., 'B' must be MxN or NxM).
if ~ismatrix(B); error('''B'' must be a vector or a matrix'); end
[bM,bN] = size(B);
if  bM ~= N && bN ~= N 
    error('The size of ''B'' must match the size of ''A'' in at least one dimension');
end

% rotate 'A' so it is 1xN and 'B' so that it is MxN, so that averaging 
% occurs along rows (for consitency later).
rotA = false;       % keep track of original orientation
rotB = false;       % keep track of original orientation
if iscolumn(A)
    A = A';
    rotA = true;    % keep track of change in orientation  
end
if bM == N
    B = B';
    rotB = true;    % keep track of change in orientation    
end
% *** NOTE: If 'B' is a square matrix, this could end up causing problems 
% if averaging is meant to occur verticall.  Additional support for such a 
% use-case should be implemented here.


% if data_range is not specified, define it
if nargin < 4 || isempty(data_range); data_range = [min(A), max(A)]; end

% Check that data_range is a two-element vector
if length(data_range) ~= 2
    error('''data_range'' must be a 2-element vector'); 
end
data_range = sort(data_range);


% if func isn't specified, define it
if nargin < 5; func = @nanmean; end

%% Define Bins

% Bin size and edges
bin_edges = data_range(1):binsize:data_range(end);

% For each data in 'a', find the corresponding bin number
bins = discretize(A, bin_edges)';
% Identify NaN bins (which occur outside of data_range of interest)
nan_bins = isnan(bins);
% Reduce dimenion of 'a' and 'b' as appropriate
bins = bins(~nan_bins);
A = A(~nan_bins);
B = B(:,~nan_bins);
[M,N] = size(B);  % overwrite previous value for N

%% Bin average (ignoring NaN's)

% for 'A' vector:
avgA = accumarray(bins,A,[],func,NaN);

% for 'B': first check if B is a vector or a matrix
if isvector(B)
    avgB = accumarray(bins,B,[],func,NaN);
elseif ismatrix(B)
    [xx,yy] = ndgrid(1:M, bins );
    subs = [xx(:) yy(:)];
    val = reshape( B, [1,M*N] );
    avgB = accumarray(subs,val,[],func,NaN);
end

%% Clean up data

% rotate 'avgA' and 'avgB' to match the original orientations of 'A' and 'B'
if rotA; avgA = avgA'; end
if rotB; avgB = avgB'; end