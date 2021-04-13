function [scHand,ctRangeSpecFlag] = TS_plot(CT,SA,ColVar,CT_range,SA_range,Col_range)
% TS_PLOT creates a temperature-salitiny scatter plot
%
%   TS_plot(CT,SA) creates a scatter plot in temperature-salitiny space
%   with background contours of sigma (potential density anomaly). CT is
%   conservative temperature and SA is absolute salinity (as per TEOS10). CT
%   and SA must be arrays of the same size (matrices are permitted).
%
%   TS_plot(CT,SA,ColVar) uses a third variable (e.g. pressure) to colour
%   the the points. ColVar must be of the same size as CT and SA, or a
%   vector with one side length matching CT and SA (e.g, CT,SA may be
%   matrices in time and depth, and ColVar could be a time vector).
%
%   TS_plot(CT,SA,ColVar,CT_range,SA_range,Col_range) optionally defines
%   the plot extents in each of the three variables.
%
%   Use empty inputs to use default values. For example, create a T-S plot
%   with a specified SA_range, but no colour-defining variable and default
%   CT_range: TS_plot(CT,SA,[],[],SA_range);
%
%   scHand = TS_plot(...) returns the plot handle for the scatterplot.
%
%   This functions requires the TEOS10 toolbox to calculate sigma contours.
%   It is available for free at: http://www.teos-10.org/software.htm
%
%   S.D.Brenner, 2018
%
%   Updated:
%       Sept., 2019: Added output argument to get scatter plot handle
%       Oct., 2019: Added ability for ColVar to be a vector with one dimension
%           matching CT,SA
%       Oct., 2019: Added freezing temperature line
%       Oct., 2019: Fixed default CT_range lower bound to account for negative
%          temperature values (maybe?)

%% Check number of inputs and fill unused inputs with default values

% Throw an error if CT and SA are not both specified
if nargin < 2
    error('Must specify at least CT and SA');
end

% If ColVar is not specified, replace with empty (useful later)
if nargin < 3; ColVar = []; end
    
    
% Default values for plotting ranges
CT_range_default = roundRangeTenths( CT(:) );
SA_range_default = roundRangeTenths( SA(:) );
Col_range_default = roundRangeTenths( ColVar(:) );

% Replace variables by default values if they're not specified
colRangeSpecFlag = false;
ctRangeSpecFlag = false;
saRangeSpecFlag = false;
if nargin < 6 || isempty(Col_range); Col_range = Col_range_default;
else colRangeSpecFlag = true; end
if nargin < 5 || isempty(SA_range); SA_range  = SA_range_default;
else saRangeSpecFlag = true; end
if nargin < 4 || isempty(CT_range); CT_range  = CT_range_default;
else ctRangeSpecFlag = true; end



%% Organize data and error checking

% Check that CT and SA are the same size and throw error if they aren't
if any( size(CT)~=size(SA) )
    error('CT and SA must be of the same size');
end

% If SA and CT are arrays, and if ColVar is specified, then it should
% either be the same size as SA and CT in at least one dimension.
% Then if it is only the same in one dimension, it will need to be repeated
% to be the same size in total.
% This piece of checking and repetition is not currently robust
if size(CT,2) == size(ColVar,2) && size(CT,1) ~= size(ColVar,1)
    ColVar = repmat( ColVar, [ size(CT,1), 1 ] );
elseif size(CT,1) == size(ColVar,1) && size(CT,2) ~= size(ColVar,2)
    ColVar = repmat( ColVar, [1, size(CT,2) ] );  
end  
% If ColVar is specified check that it is also the same size
if any( size(ColVar)~=size(CT) ) && ~isempty(ColVar)
    warning('Size of ColVar variable does not match size of CT, SA.  ColVar not used');
end
    
% Change data varibles into column vectors (accounts for matrix inputs)
CT = CT(:);
SA = SA(:);
ColVar = ColVar(:);    

% Sort data based on ColVar if it exists and is of the correct size
if ~isempty(ColVar) && length(ColVar) == length(CT)
    [ColVar,idx] = sort(ColVar);
    CT = CT(idx);
    SA = SA(idx);
end
 
% If ColVar is not specified or the wrong length, replace with grey
if  isempty(ColVar) || length(ColVar)~=length(CT)
    ColVar = [0.5,0.5,0.5];
    Col_range = [0,1];
end


%% Modify CT_range lower bound
% If CT_range is not specified, but the freezing line is being included,
% extend the lower bound of the CT_range so that the line is visible
% (this extension must happen here in the code so that the sigma curves
% cover the whole domain, even though the actually line is added at the
% end)

sa = linspace( SA_range(1),SA_range(2) );
ctFr = gsw_CT_freezing(sa,0);

if ~ctRangeSpecFlag
    ctFrRound = roundRangeTenths(ctFr);
    CT_range(1) = ctFrRound(1) ;
end

% If CT_range IS specified, then ctFr must be contained within it:
if ctRangeSpecFlag
    ctFr(ctFr <= min(CT_range)) = NaN;
    ctFr(ctFr >= max(CT_range)) = NaN;
end



%% Generate curves of CT vs SA:
grey = 0.7*[1,1,1];


ct = linspace( CT_range(1),CT_range(2) );
[sa_mat,ct_mat] = meshgrid(sa,ct);
sigma = gsw_sigma0(sa_mat,ct_mat);

[c,hContour] = contour(sa,ct,sigma,'color',grey,'Levelstep',0.5); 
clabel(c,hContour,'color',grey,'interpreter','latex') % default labels
addlistener(hContour,'MarkedClean',@(a,b)ReFormatText(a));


%% Overlay T-S Points onto sigma curves

sz = 10; % point size
hold on;
scHand =  scatter(SA,CT,sz,ColVar,....
        'filled',...
        'MarkerEdgeColor','none');
hold off;    
set(gca,'xlim',SA_range,...
        'ylim',CT_range,...
        'clim',Col_range);
        
% Labelling
xlabel('S_A [g/kg]');
ylabel(['\Theta [',char(176),'C]']);



%% Add CT_freezing line 


hold on;
plot( sa, ctFr, 'k--', 'linewidth',1);
hold off;

end


function roundX = roundRangeTenths(X)
    % should be robust to negative numbers
    minX = min(X(:));
    maxX = max(X(:));
    minRoundX = ( floor(10*minX ) - floor( minX )/10 )/10;
    maxRoundX = (  ceil(10*maxX ) +  ceil( maxX )/10 )/10;
    roundX = [ minRoundX, maxRoundX ];
end




function ReFormatText(h)
  % get all the text items from the contour
  t = get(h,'TextPrims');
  for ii=1:length(t)
    % get the current value (Matlab changes this back when it 
    %   redraws the plot)
%     v = str2double(get(t(ii),'String'));
    str = get(t(ii),'String');
    % Update with the format that contains sigma_theta
    set(t(ii),'String', strcat(str,'$\sigma_{_\Theta}$') );
  end
end  