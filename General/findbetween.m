function idx = findbetween(x,xRange)
% FINDBETWEEN finds indices of values contained between two bounds
%
%   idx = findbetween(x,xRange) returns the indices in x that are within
%   the range of xRange (inclusive of endpoints)
%
%   S.D.Brenner, 2019
     
%% Error checking

if length(xRange) ~=2 || ~issorted(xRange) || ~isnumeric(xRange) || any(isnan(xRange))
    error('Value xRange must be a 1x2 vector of numeric type in which the second element is larger than the first and may be Inf');
end
if ~isvector(x)
    error('x must be a vector');
end

%% Actual function

idx = find( x >= xRange(1) & x <= xRange(2) );
    
end
    
    
    