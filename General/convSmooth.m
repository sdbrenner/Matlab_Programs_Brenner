function [Zs] = convSmooth(Z,N,sigma)
% [Zs] = CONVSMOOTH(x,y,Z)
%
%   2-dimensional convolution smoother using a N-point guassian kernel
%   See:    https://www.mathworks.com/help/matlab/data_analysis/convolution-filter-to-smooth-data.html
%   and:    https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
%   This may also be a poor-man's version of imgaussfilt ?


% Set defaults
sigma = 1;
N = 9; % (there is supposed to be some relationship between sigma and N)


% Build N-point Gaussian smoothing Kernal:
x = linspace( -3*sigma, 3*sigma, N);
[X,Y] = meshgrid( x , x );
K = exp( -(X.^2 + Y.^2) / sqrt(2) );
K = exp( -(X.^2 + Y.^2) / (2*sigma^2) );
K = K./sum(K(:)); % normalize weights to sum to 1

% To reduce boundary effects, mirror the input matrix across all bondaries
Zcent = [fliplr(Z), Z, fliplr(Z)];
Ztopbot = flipud(Zcent);
Zfull = [Ztopbot; Zcent; Ztopbot];

% Apply convolution
Zconv = conv2(Zfull,K,'same');

% Extract original matrix subset
[numRow,numCol] = size(Z);
extractRows = (1:numRow) + numRow;
extractCols = (1:numCol) + numCol;
Zs = Zconv( extractRows, extractCols);

end

