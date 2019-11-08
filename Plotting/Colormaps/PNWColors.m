function cmap = PNWColors(Name,M)
% PNWColors Colormaps inspired by the Pacific Northwest. 
%
%   cmap = PNWColors(Name,M) generates an M-by-3 matrix, cmap, of rbg
%   colormap data corresponding to a given colormap Name, from the list of
%   PNW Colors produced by Jake Lawlor.  If M is not specified, the default
%   colormap length will be used (5 to 8).  Note that the colormap name is
%   case sensitive.
%
%   PNWColors called without any inputs will generate a new plot figure
%   showing all of the colormap options and names.
%   
%   Colormaps were copied from:
%      https://github.com/jakelawlor/PNWColors
%   This code simply converts his beautiful work to something useable in 
%   Matlab.
%
%   Example usage: 
%       surf( peaks );
%       shading interp;
%       cmap = PNWColors( 'Sailboat', 32 );
%       colormap( cmap );
%       colorbar;
%
%   S.D.Brenner, 2019

%% First, build the colormaps, then check for input requests


% Define colormaps in hex code
% ( these hex-keys are copy-pasted directly from
%   https://github.com/jakelawlor/PNWColors/blob/master/R/PNWColors.R )
hexCol.Starfish = {'#24492e', '#015b58', '#2c6184', '#59629b', '#89689d', '#ba7999', '#e69b99'};
hexCol.Shuksan = {'#33271e', '#74677e', '#ac8eab', '#d7b1c5', '#ebbdc8', '#f2cec7', '#f8e3d1', '#fefbe9'};
hexCol.Bay = {'#00496f', '#0f85a0', '#edd746', '#ed8b00', '#dd4124'};
hexCol.Winter = {'#2d2926', '#33454e', '#537380', '#81a9ad', '#ececec'};
hexCol.Lake = {'#362904', '#54450f', '#45681e', '#4a9152', '#64a8a8', '#85b6ce', '#cde5f9', '#eef3ff'};
hexCol.Sunset = {'#41476b', '#675478', '#9e6374', '#c67b6f', '#de9b71', '#efbc82', '#fbdfa2'};
hexCol.Shuksan2 = {'#5d74a5', '#b0cbe7', '#fef7c7', '#eba07e', '#a8554e'};
hexCol.Cascades = {'#2d4030','#516823','#dec000','#e2e260','#677e8e','#88a2b9'};
hexCol.Sailboat = {'#6e7cb9', '#7bbcd5', '#d0e2af', '#f5db99', '#e89c81', '#d2848d'};
hexCol.Moth = {'#4a3a3b', '#984136', '#c26a7a', '#ecc0a1', '#f0f0e4'};

% Convert to Matlab compatible colormaps
flds = fields(hexCol);
N = length(flds);
for n = 1:N
    fldName = flds{n};
    cmap.(fldName) = hex2col( hexCol.(fldName) );
end

%% Parse inputs


% If no arguements are specified, show all color options
if nargin == 0
    makePlots(cmap); 
    clear cmap
    return;
else
    % Check for existing name:
    if ~isfield(cmap,Name)
        makePlots(cmap); 
        error('Pick an existing colormaps');
    end
    cmap = cmap.(Name);
end

% If a specific number of colours are specified, adjust the colormap:
if nargin > 1 && isnumeric(M)
    cmap = cinterp2( cmap, M );
end


end

%% EMBEDDED FUNCTIONS %% ==================================================

function col = hex2col(S)
    
    for n = 1:length(S)
        s = S{n};
        col(n,1) = hex2dec( s(2:3) )/255;
        col(n,2) = hex2dec( s(4:5) )/255;
        col(n,3) = hex2dec( s(6:7) )/255;
    end
end
    

function makePlots(cmap)

    flds = fields(cmap);
    N = length(flds);
    figure; clf;
    for n = 1:N
        fldName = flds{n};
        plotCbar(cmap.(fldName),fldName,n,N)
    end    
    set(gcf,'color','w');
end



function plotCbar(col,name,n,N)
    ypos = 0.1;
    height = 0.8;
    xpos = n/(N+2);
    width = 0.7/(N+2);

    ax = subplot(1,N,n);
    colormap(ax,col);
    axis off
    cb = colorbar;
    
    cb.Position = [ xpos, ypos, width, height ];
    cb.Ticks = [];
    title(cb,name);
end



function cmapM = cinterp2(cmap,M)
    
    L = length(cmap);

    x = 1:3;
    y = 1:L;
    [X,Y] = meshgrid(x,y);
    
    yq = linspace(1,L,M);
    [Xq,Yq] = meshgrid(x,yq);
    
    cmapM = interp2(X,Y,cmap,Xq,Yq);

end
