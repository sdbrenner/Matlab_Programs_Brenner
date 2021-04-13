function colour = cinterp(data, clim, cmap)
% CINTERP colormap interpolation
% 
%	colour = CINTERP(data, clim, cmap)
% 
%	Colour interpolation: Given a colormap (cmap) with limits defined by
%	clim (i.e: colormap(cmap); caxis(clim) ), cinterp will linearly
%	interpolate cmap to find the colours that correspond to data values
%	within clim
%
%   Usage example: 
%   colour lines in a plot based on an additional variable, cData
%       cmap = jet(20);
%   	clim = [0,10];
%       cData = [1,2,3,7];
%       cols =  cinterp( cData, clim, cmap );
%       figure;
%       hold on;
%       set(gca,'ColorOrder',cols);
%       plot( magic(4) );
%       caxis(clim);
%       colormap(cmap);
%       colorbar;
%
%	S.D.Brenner, 2018


% Number of distinct colours in map:
N = size(cmap,1);
linear_color_vect = linspace(clim(1),clim(2),N);

% interpolate each vector of colour individually
c1 = interp1( linear_color_vect, cmap(:,1), data );
c2 = interp1( linear_color_vect, cmap(:,2), data );
c3 = interp1( linear_color_vect, cmap(:,3), data );

% the final colour is a vector of the three colours
colour = [c1(:),c2(:),c3(:)];

end