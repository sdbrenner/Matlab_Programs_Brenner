function [dux,duy] = m_gradient(u,lon,lat)
% M_GRADIENT approximates spatial gradients on a lat-lon grid

    % Transform from lat, lon to distance in [m]:
%     earthRad = 6378.137;                       % [km] (equitorial radius) 
    earthRad = 6371.009;                        % [km] (average radius) 
    deg2meter = deg2rad(1) * earthRad * 1e3;    % [meters per degree latitude]
    lon0 = median(lon(:));
    y = deg2meter*(lat);                        % [m]
    x = deg2meter*cosd(lat).*(lon-lon0);        % [m]

    % Calculate grandient wrt grid i-j idx grid
    [dxi,dxj]=gradient(x);
    [dyi,dyj]=gradient(y);
    [dui,duj]=gradient(u);
    
    % Determinant of Jacobian
    DET = dxi.*dyj- dxj.*dyi;
    if abs(DET)<1e-12, error('bad transformation'); end
   
    % Transform the gradient [dui,duj] into dux, duy:
    dux = (1./DET) .* ( dui.*dyj - duj.*dyi );      
    duy = (1./DET) .* ( duj.*dxi - dui.*dxj );
    
end