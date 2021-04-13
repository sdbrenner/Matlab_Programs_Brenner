function [xBin,yBin] = binFun(X,varargin)
% BINFUN Apply function to data within bins that are defined in terms of a
% dependent variable (e.g., bin-averaging).
%
%   xBin = binFun(x)
%   
%   [xBin,yBin] = binFun(x,y)
%   [xBin,yBin] = binFun(x,y,binSize) 
%   [xBin,yBin] = binFun(x,y,binEdges) 
%   [xBin,yBin] = binFun(...,Name,Value) 
%
%   S.D.Brenner, 2020.

%% Parse arguments
[x,y,binEdges,fun,fillValue,opts] = parseInputs(X,varargin{:});


%% Discretize

bins = discretize( x, binEdges );
bins = bins(:);

% remove values where discretize has output NaN values
% (these represent values of x outside of the bounds of binEdges)
x( isnan(bins) ) = [];
if ~isempty(y)
    y( isnan(bins),:) = [];
end
bins( isnan(bins) ) = [];


%% Apply binning function

numBins = length(binEdges)-1;
% Use accumarray to perform function on 'x' values in each bin 
xBin = accumarray( bins, x, [numBins,1], fun, fillValue );



% Make application for 'y' more general so it applies for both vector and
% matrix inputs (also works fine if y is empty)
numCol = size(y,2);
[a,b] = ndgrid( 1:numCol , bins );
binsPointPairs = [a(:),b(:)];
yBin = accumarray( binsPointPairs , reshape(y.',1,[]) , [numCol,numBins],...
                   fun , fillValue ).';


%% Adjust outputs based on other options

% Check if returning calculated xBin or binCenters
binCents = mean( [binEdges(1:end-1);binEdges(2:end) ] );
if strcmp(opts.returnX,'binCenters')
    xBin = binCents(:);
end

% Check if missing values should be returned or excluded
existFun = @(x) ~any(isempty(x));
missingInd = ~accumarray(bins, x, [numBins,1], existFun );
if ~opts.includeMissing
    xBin( missingInd ) = [];
    yBin( missingInd,:) = [];
end


% Check if x or y outputs need to be rotated (they should be returned in
% the same orientation that they were input)
if opts.rotX; xBin = xBin.'; end
if opts.rotY; yBin = yBin.'; end



end


%%
function [x,y,binEdges,fun,fillValue,opts] = parseInputs(x,varargin)

%% Set up and run input parser

p = inputParser;

% Define default values
defaultBinFun = @(x) mean(x,'omitnan');
defaultBinLimits = [ min(x) , max(x) ];

% Define expected strings
expectedBinMethods = {'auto','scott','fd','integers','sturges','sqrt'};
expectedReturnX = {'binCalc','binCenters'};

% Define validation functions
validX =         @(x) validateattributes(x,{'numeric'},{'vector'});
validY =         @(x) validateattributes(x,{'numeric'},{'2d'});
validBinInput =  @(x)validateattributes(x,{'numeric'},{'vector','increasing'});
validBinFun =    @(x) validateattributes(x,{'function_handle'},{});
validFillValue = @(x) validateattributes(x,{'numeric'},{'scalar'});
validBinLimits = @(x) validateattributes(x,{'numeric'},{'real','numel',2,'increasing'});
validNumBins =   @(x) validateattributes(x,{'numeric'},{'real','scalar','nonnegative'});
validBinMethod = @(x) any( validatestring(x,expectedBinMethods) );
validReturnX =   @(x) any( validatestring(x,expectedReturnX) );
validIncMssng =  @(x) validateattributes(x,{'logical'},{});


% Add attributes (required/optional):
addRequired(p,'x',validX);
addOptional(p,'y',[],validY);
addOptional(p,'binInput',[],validBinInput);

% Add parameters (Name-Value pairs)::
addParameter(p,'binFunction',defaultBinFun,validBinFun);
addParameter(p,'fillValue',NaN,validFillValue);
addParameter(p,'binLimits',defaultBinLimits,validBinLimits);
addParameter(p,'numBins',[],validNumBins);
addParameter(p,'binMethod','auto',validBinMethod);
addParameter(p,'returnX','binCalc',validReturnX);
addParameter(p,'includeMissing',true,validIncMssng); 
                                 
% Run input parser  
parse(p,x,varargin{:});                            
       
   

%% Extract variables from parser

x = p.Results.x;
y = p.Results.y;
binInput = p.Results.binInput;
fun = p.Results.binFunction;
fillValue = p.Results.fillValue;

binLimits = p.Results.binLimits;
numBins = p.Results.numBins;
binMethod = p.Results.binMethod;

opts.returnX = p.Results.returnX;
opts.includeMissing = p.Results.includeMissing;

%% Find or generate vector of bin edges

% If no input is given for binInput, generate bin edges using either
% linspace or histcounts:
if isempty(binInput)
    if isdefault(p,'binLimits') && isdefault(p,'numBins')
        [~,binEdges] = histcounts( x,'BinMethod', binMethod );
    elseif isdefault(p,'binLimits') && ~isdefault(p,'numBins')
        [~,binEdges] = histcounts( x, numBins,'BinMethod', binMethod );
    elseif ~isdefault(p,'binLimits') && isdefault(p,'numBins')
        [~,binEdges] = histcounts( x,'BinMethod', binMethod,...
                                     'BinLimits',binLimits);
    else
        binEdges = linspace( binLimits(1), binLimits(2), numBins+1 );
    end
    
% If binInput is a scalar, take it as the bin size:    
elseif isscalar(binInput)
    binSize = binInput;
    binEdges = unique([binLimits(1):binSize:binLimits(2),binLimits(2)]);
    
% If binInput is a vector, take it as the bin edges directly:    
elseif isvector(binInput)
    binEdges = binInput;
% All options should now be accounted for, so I can't see why this error
% would ever appear:
else
    error('Unexpected binInput error');
end


%% Rotate vectors (if necessary)

% Rotate 'x' so it is a column vector with size Nx1
opts.rotX = isrow(x);
x = x(:);

% If necessary, rotate 'y' so it is either:
%    (1) a column vector with size Nx1
%    (2) a matrix with size NxM
opts.rotY = 0;
if ~isempty(y)
    if size(y,1) == length(x)

    elseif size(y,2) == length(x)
        opts.rotY = 1;
        y = y.';
    else
        error('The size of ''y'' must match the size of ''x'' in at least one dimension');    
    end
end

% Rotate binEdges so it is a row vector
binEdges = binEdges(:).';

end
%%


function value = isdefault(p,Name)
    value = any(  strcmp(p.UsingDefaults, Name)  );
end