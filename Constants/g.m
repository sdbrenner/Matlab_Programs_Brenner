function grav_const = g(L,T)
% G returns the values of the standard values used for  acceleration due to 
% gravity on Earth
%
%   g = 9.81 m/s/s
%
%   g(L,T) can return g in other length/time units consistent with a
%   'symunit' collection.  L are the units associated with length (default
%   'm'), and T are the units associated with time (default 's');
%
%   S.D.Brenner, 2019

% Default units:
su = symunit;
if nargin < 2 || isempty(T); T = 's'; end
if nargin < 1 || isempty(L); L = 'm'; end

% Unit conversion factor (if other units are requested)
unitConv = unitConversionFactor( su.m/su.s/su.s,...
                                  su.(L)/su.(T)/su.(T)   );

% Value                                          
grav_const = 9.81 * double(unitConv);


end