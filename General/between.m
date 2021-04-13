function idx = between(x,xRange)
% BETWEEN finds indices of values contained between two bounds
%
%   idx = between(x,xRange) returns the logical 1 (TRUE) for in x that are
%   within the range of xRange and logical 0 (FALSE) for values outside of
%   the range (inclusive of endpoints)
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

idx = ( x >= xRange(1) & x <= xRange(2) );
    
end
    
    
    