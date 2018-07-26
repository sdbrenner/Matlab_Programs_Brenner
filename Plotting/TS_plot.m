function TS_plot(CT,SA,ColVar,CT_range,SA_range,Col_range)
%
% TS_plot(CT,SA) creates a scatter plot in temperature-salitiny space with
% background contours of sigma (potential density anomaly)
% CT is conservative temperature and SA is absolute salinity (as per
% TEOS10)
% CT and SA must be arrays of the same size (matrices are permitted).  
%
% TS_plot(CT,SA,ColVar) uses a third variable (e.g. pressure) to colour the
% the points.  
% ColVar must be of the same size as CT and SA.
%
% Plot extents can be optionally defined using 
% TS_plot(CT,SA,ColVar,CT_range,SA_range,Col_range)
%
% Use empty inputs to use default values. 
% For example, create a T-S plot with a specified SA_range, but no 
% colour-defining variable and default CT_range:
% TS_plot(CT,SA,[],[],SA_range);
%
% This functions requires the TEOS10 toolbox to calculate sigma contours.
% It is available for free at: http://www.teos-10.org/software.htm
%
% Samuel Brenner, 2018


%% Check number of inputs and fill unused inputs with empty [];
if nargin < 2
    error('Must specify at least CT and SA');
elseif nargin == 2
    ColVar = [];
    CT_range = [];
    SA_range= [];
    Col_range = [];
elseif nargin == 3
    CT_range = [];
    SA_range= [];
    Col_range = [];
elseif nargin == 4
    SA_range= [];
    Col_range = [];
elseif nargin == 5
    Col_range = []; 
end

%% Organize data and error checking

% Check that CT and SA are the same size and throw error if they don't
% match
if any( size(CT)~=size(SA) )
    error('CT and SA must be of the same size');
end

% If ColVar is non empty, check that it is also the same size.
if any( size(ColVar)~=size(CT) ) && ~isempty(ColVar)
    warning('Size of ColVar variable does not match size of CT, SA.  ColVar not used');
end
    
% Change data varibles into column vectors (accounts for matrix inputs).
CT = CT(:);
SA = SA(:);
ColVar = ColVar(:);    

% Sort data based on ColVar if it exists and is of the correct size
if ~isempty(ColVar) && length(ColVar) == length(CT)
    [ColVar,idx] = sort(ColVar);
    CT = CT(idx);
    SA = SA(idx);
end


%% Set default values:

% If CT_range, SA_range, and Col_range aren't specified, set the range to 
% go from 90% of the lowest data point to 110% of the highest data point, 
% rounded up or down as appropriate to the nearest tenth.
if isempty(CT_range)
    CT_range = [floor(9*min(CT))/10 ceil(11*max(CT))/10 ];
end
if isempty(SA_range)
    SA_range = [floor(9*min(SA))/10 ceil(11*max(SA))/10 ];
end   
if isempty(Col_range) && ~isempty(ColVar)
    Col_range = [floor(9*min(ColVar))/10 ceil(11*max(ColVar))/10 ];
elseif isempty(Col_range) && isempty(ColVar)
    Col_range = [0,1];
end    

if isempty(ColVar) || length(ColVar)~=length(CT)
    % if ColVar is not specified or of the wrong length, replace with
    % 'black'
    ColVar = [0,0,0];
end

%% Generate curves of CT vs SA:

sa = linspace( SA_range(1),SA_range(2) );
ct = linspace( CT_range(1),CT_range(2) );
[sa_mat,ct_mat] = meshgrid(sa,ct);
sigma = gsw_sigma0(sa_mat,ct_mat);

[c,hContour] = contour(sa,ct,sigma,'k'); 
clabel(c,hContour) % default labels



%% Overlay T-S Points onto sigma curves

sz = 5; % point size.
hold on;
scatter(SA,CT,sz,ColVar,....
        'filled',...
        'MarkerEdgeColor','none');
set(gca,'xlim',SA_range,...
        'ylim',CT_range,...
        'clim',Col_range);
        
%% Labelling

xlabel('S_A');
ylabel('\Theta');

% modify hidden properties of curves to include '?' symbol
% (apparently this doesn't work inside the function?)
N = length( hContour.TextPrims );
for n = 1:N
    hContour.TextPrims(n).String = [hContour.TextPrims(n).String,'?'];
end
hold off;

end
