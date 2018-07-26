function [U,V,err] = READ_icemotion(filename)
% Read in icemotion files from the Polar Pathfinder data product
% https://nsidc.org/data/nsidc-0116
% Output velocities U and V, and error std are in units of [km/day]

fileID = fopen(filename);
[~,~,machinefmt,~] = fopen(fileID);

L = 361^2;
A = fread(fileID,'int16',machinefmt);
B = reshape(A,[3,L]);


% [cm/sec] ? Scaled by a factor of 10
U = reshape(B(1,:),[361,361])/10;
V = -reshape(B(2,:),[361,361])/10;

% The third variable contains the square root of the estimated error 
% variance, scaled by a factor of 10, at a given location:
err = reshape(B(3,:),[361,361])/10;

% Convert output variables to units of [km/day] instead of [cm/s]:
U = U * 86400/100000;
V = V * 86400/100000;
err = err * 86400/100000;

end