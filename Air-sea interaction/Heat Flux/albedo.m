function alpha = albedo(mu)
% alpha = ALBEDO(mu) finds the ocean surface albedo parameterized as a 
%   function of solar zenith angle.
%
%   'mu' is the cosine of the local solar zenith angle.
%
%   from: 
%   "Comparison of Regional Clear-Sky Albedos Inferred from Satellite  
%   Observations and Model Computations",
%   B. P. Briegleb et. al., 1986
%   doi:  10.1175/1520-0450(1986)025<0214:CORCSA>2.0.CO;2
%
%   Samuel Brenner, 2018



alpha = 2.6./(mu.^1.7 + 0.065) + 15.0*(mu-1).*(mu-0.5).*(mu-1.0);
alpha = real(alpha)/100;
alpha(mu<0) = NaN; % albedo is not meaningful during nighttime when there is no solar radiation

end