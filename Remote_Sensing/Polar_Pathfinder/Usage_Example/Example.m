%% Usage example
% This usage example require the Matlab package 'm_map', produced by Rich
% Pawlowicz and available for free at https://www.eoas.ubc.ca/~rich/map.html


%% READ_icemotion.m
% Read Polar Pathfinder icemotion binary file and plot icemotion vectors

% Build Polar Pathfinder coordinates
% * note: a text file containing EASEgrid coordinates is available at 
% https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0116_icemotion_vectors_v3/tools/
load('icemotion_coords.mat');


% Idenitify file and read in using 'READ_icemotion'
filename = 'icemotion.grid.week.2014.39.n.v3.bin';
[U,V,err] = READ_icemotion(filename);

% Plot vector directions
figure(1);
quiver(X,Y,U,V);
set(gca,'xlim',[90,270],'ylim',[90,270]);
axis square;


%% map2EASEgrid.m
% Load coastline file and convert to EASEgrid coordinates

% Load file
load coast;

% Convert lat,lon map coordinates to EASEgird coordinates
[coast_x,coast_y] = map2EASEgrid(long,lat);

% Add to coastlines to current plot
hold on;
plot(coast_x,coast_y,'k');


%% EASEgrid2map.m
% Show icemotion speed on a map

% Convert U,V cartesian vectors to scalar speed value
Spd = sqrt( U.^2 + V.^2 );

% Convert EASEgrid X,Y coordinates to lat, lon map coordinates
[EASE_lon,EASE_lat] = EASEgrid2map(X,Y);
% * note: the coordinates file includes lat,lon coordinate values.  This 
%   conversion is only done as an example, and isn't strictly necessary


% Create figure with mapping projection
figure(2);
m_proj('stereographic','lat',90,'long',30,'radius',25);

% Plot ice speed
m_pcolor(EASE_lon,EASE_lat,Spd);

% Add coastlines and grid
m_coast('patch',[0.8,0.8,0.8],'edgecolor','none');
m_grid;


%% vec_EASEgrid2map.m
% Show icemotion vectors on a map

% Convert icemotion vectors from EASEgrid to mapping coordinates
[EASE_lon,EASE_lat,U_map,V_map] = vec_EASEgrid2map(X,Y,U,V);

% Create figure with mapping projection
figure(3);
m_proj('stereographic','lat',90,'long',30,'radius',25);

% Plot ice speed
m_quiver(EASE_lon,EASE_lat,U_map,V_map);

% Add coastlines and grid
m_coast('patch',[0.8,0.8,0.8],'edgecolor','none');
m_grid;




