function grav_const = g(L,T)
% G returns the values of the standard values used for  acceleration due to 
% gravity on Earth
%
%   g = 9.81 m/s/s
%
%   g(L,T) can return g in other units

su = symunit;
if nargin <2; L = 'm'; T = 's'; end


unitConv = unitConversionFactor( su.m/su.s/su.s,...
                                          su.(L)/su.(T)/su.(T)   );                                     
grav_const = 9.81 * double(unitConv);


end