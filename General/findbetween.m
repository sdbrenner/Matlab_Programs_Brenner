function idx = findbetween(x,x_range)
    % idx = between(x,x_range) returns the indices in x that are within the
    % range of x_range (nclusive of endpoints)
    
    
    % Error checking
    if length(x_range) ~=2 
        error('x_range must be length 2');
    end
    if ~isvector(x)
       error('x must be a vector');
    end
    
    % Actual function
    
    idx = find( x >= x_range(1) & x <= x_range(2) );
    
end
    
    
    