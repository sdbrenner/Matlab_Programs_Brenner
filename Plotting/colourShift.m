function c  = colourShift(col,scale,direction)
% COLOURSHIFT shifts a colormap lighter or darker
%
%   c  = colourShift(col,scale,direction) where 'col' is the colormap to
%   shift, 'scale' is the amount of shift (scalar value between 0 and 1),
%   and 'direction' is one of either 'light' or 'dark', specifying the
%   direction of the shift
%   
%   S.D.Brenner, 2021


if isempty(scale)
    scale = 0.85;
end


switch direction
    case 'light'
        c = 1-scale*(1-col);
    case 'dark'
        c = scale*col;
    otherwise
        error("'direction' must be one of either 'light' or 'dark'")
end

end