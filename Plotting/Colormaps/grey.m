function col = grey( intensity )
% GREY colour maps
%
%   col = grey( intensity ) returns an Nx3 matrix of colormap values, for N
%   the length of the intensity vectro


    col = intensity(:)*[1,1,1];
end